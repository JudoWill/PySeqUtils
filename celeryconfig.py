from kombu import Exchange, Queue
BROKER_URL = 'amqp://'
CELERY_RESULT_BACKEND = 'amqp://'

CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'

CELERY_IMPORTS = ('HIVTransTool',
                  'TreeingTools',)

CELERY_QUEUES = [Queue('HIVTransTool'),
                 Queue('default'),
                 Queue('long-running')]

CELERY_ANNOTATIONS = {
   'HIVTransTool.map_seqs_to_ref': {'rate_limit': '10/m'}
}