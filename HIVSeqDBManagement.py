__author__ = 'will'
from ming import create_datastore
from ming import Session
from ming import Field, schema
from ming.declarative import Document
from GeneralSeqTools import fasta_reader
from fileinput import FileInput
from HIVTransTool import map_seqs_to_ref, load_region, region_linker
from copy import deepcopy
from itertools import islice

from celery import Celery
from celery import group
import logging

celery = Celery()
celery.config_from_object('celeryconfig')

bind = create_datastore('testsequence')
session = Session(bind)


class SequenceObject(Document):

    class __mongometa__:
        session = session
        name = 'sequence_object'

    _id = Field(schema.ObjectId)
    PatID = Field(str)
    NucSeq = Field(str)
    AASeq = Field(str)
    RegionName = Field(str)
    IsRaw = Field(bool)
    BeenTranslated = Field(bool)
    Source = Field(str)


def load_raw_seqs_to_db(seq_iterable, source, is_nuc=True, RegionName='Genome'):

    field = 'NucSeq' if is_nuc else 'AASeq'
    for num, (name, seq) in enumerate(seq_iterable):
        if (num == 10) | (num % 10000 == 0):
            print 'Integrated: %i sequences!' % num
        obj = SequenceObject({
            'PatID': name,
            field: seq,
            'IsRaw': True,
            'Source': source,
            'RegionName': RegionName,
            'BeenTranslated': False
        })
        obj.m.save()


def load_fasta_to_db(fasta_files, source, is_nuc=True, RegionName='Genome'):

    seq_iterable = fasta_reader(FileInput(fasta_files))
    load_raw_seqs_to_db(seq_iterable, source, is_nuc=is_nuc, RegionName=RegionName)


def get_seqs_for_LANL_translation():

    needed = SequenceObject.m.find({'BeenTranslated': False}).count()
    print 'Need to translate %i sequences' % needed
    print 'doing query'
    query = SequenceObject.m.find({'BeenTranslated': False}).all()

    print 'yielding seqs'
    for row in query:
        yield (row['PatID'], row['NucSeq'])


@celery.task(name='HIVSeqDBManagement.query_LANL',
             queue='HIVSeqDBManagement')
def query_LANL(input_seqs):

    region_dict = load_region()
    update_names = set()

    for row in map_seqs_to_ref(input_seqs):
        update_names.add(row['Name'])
        obj = SequenceObject({
            'PatID': row['Name'],
            'RegionName': row['RegionName'],
            'NucSeq': row['QueryNuc'],
            'AASeq': row['QueryAA'],
            'IsRaw': False,
            'BeenTranslated': True
        })
        obj.m.save()
        for region_row in region_dict[row['RegionName']]:
            nrow = region_linker(deepcopy(row), region_row)
            if nrow:
                obj = SequenceObject({
                    'PatID': nrow['Name'],
                    'RegionName': nrow['RegionName'],
                    'NucSeq': nrow['QueryNuc'],
                    'AASeq': nrow['QueryAA'],
                    'IsRaw': False,
                    'BeenTranslated': True
                })
                obj.m.save()

    logging.warning('Updating found sequences')
    for name in update_names:
        fix_doc = SequenceObject.m.find({
            'PatID': name,
            'IsRaw': True,
            'BeenTranslated': False
        }).first()
        if fix_doc:
            fix_doc['BeenTranslated'] = True
            fix_doc.m.save()

    logging.warning('DONE!')
    return len(update_names)


def take(iterable, chunk):
    return list(islice(iterable, chunk))


def yield_chunks(iterable, chunksize):
    """Yields chunks of items from an iterable."""

    chunk = take(iterable, chunksize)
    while chunk:
        yield chunk
        chunk = take(iterable, chunksize)


def annotate_database():

    wanted_seqs = get_seqs_for_LANL_translation()
    chunks = yield_chunks(wanted_seqs, 100)
    jobs = [query_LANL.subtask((chunk,)) for chunk in chunks]
    print 'generating %i jobs' % len(jobs)
    job = group(jobs)
    res = job.skew(start=0, stop=50*len(jobs)).apply_async()
    print 'waiting for results!'
    tot = 0
    for res in res.iterate():
        tot += res
        print 'Committed %i seqs' % tot


if __name__ == '__main__':

    annotate_database()