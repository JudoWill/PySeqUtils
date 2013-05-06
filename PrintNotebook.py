__author__ = 'will'
from subprocess import check_call
import shlex
import shutil
import os.path
import os
import contextlib
import argparse

from TreeingTools import tmp_directory

@contextlib.contextmanager
def tmp_changedir(new_dir):

    cur_dir = os.curdir
    try:
        os.chdir(new_dir)
        yield new_dir
    finally:
        os.chdir(cur_dir)


def make_pdf(input_notebook_path, output_pdf):
    """Converts a ipyp notebook to a pdf"""

    notebook_fname = input_notebook_path.rsplit(os.sep, 1)[-1]
    notebook_base = notebook_fname.rsplit('.', 1)[0]
    with tmp_directory() as tmp_path:
        with tmp_changedir(tmp_path):

            new_notebook_path = os.path.join(tmp_path, notebook_fname)
            new_tex_path = os.path.join(tmp_path, notebook_base+'.tex')
            new_pdf_path = os.path.join(tmp_path, notebook_base+'.pdf')
            shutil.copy(input_notebook_path, new_notebook_path)

            run_nbconvert(new_notebook_path)
            run_latex(new_tex_path)

            shutil.move(new_pdf_path, output_pdf)


def run_latex(tex_name):
    """Runs latex to convert a .tex into a pdf"""
    cmd = 'texi2pdf ' + tex_name
    cmd_list = shlex.split(cmd)
    check_call(cmd_list)


def run_nbconvert(notebook_path, convert_path = '/home/will/nbconvert'):
    """Runs the nbconvert to generate a .tex file and associated figures."""

    cmd = 'python %(nb_path)s --format=latex %(notebook)s'
    tdict = {
        'nb_path': os.path.join(convert_path, 'nbconvert.py'),
        'notebook': notebook_path
    }

    cmd_list = shlex.split(cmd % tdict)
    check_call(cmd_list)


def find_notebook(input_string):
    """Tries to interpret the input string into a path."""

    if os.path.exists(input_string):
        return os.path.abspath(input_string)
    if os.path.exists('/home/will/IpythonNotebook/' + input_string):
        return '/home/will/IpythonNotebook/' + input_string

    try:
        return find_notebook(input_string + '.ipynb')
    except IOError:
        raise IOError, 'Could not find notebook!', input_string


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A small helper program to convert notebooks to pdfs directly')
    parser.add_argument('InNotebook',
                        help='The notebook to convert.')
    parser.add_argument('-o', '--output', required=False,
                        help='The output path')
    parser.add_argument('-d', '--dropsync', action='store_true', default=False,
                        help='Copy output directly to Dropsync')

    args = parser.parse_args()

    in_path = find_notebook(args.InNotebook)
    notebook_base = in_path.rsplit(os.sep, 1)[1].rsplit('.', 1)[0]
    base_path = in_path.rsplit(os.sep, 1)[0]
    if args.output:
        out_path = args.output
    elif args.dropsync:
        out_path = '/home/will/Dropbox/Dropsync/IpythonNotebooks/' + notebook_base + '.pdf'
    else:
        out_path = in_path.replace('.ipynb', '.pdf')

    make_pdf(in_path, os.path.abspath(out_path))
