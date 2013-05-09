BROKER_URL = 'amqp://'
CELERY_RESULT_BACKEND = 'amqp://'

CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'

CELERY_IMPORTS = ('HIVTransTool',
                  'TreeingTools',)

CELERY_ANNOTATIONS = {
    'HIVTransTool.map_seqs_to_ref': {'rate_limit': '10/m'}
}