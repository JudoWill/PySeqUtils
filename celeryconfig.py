from kombu import Exchange, Queue
BROKER_URL = 'amqp://'
CELERY_RESULT_BACKEND = 'amqp://'

CELERY_TASK_SERIALIZER = 'pickle'
CELERY_RESULT_SERIALIZER = 'pickle'

BROKER_HEARTBEAT = 10

CELERY_IMPORTS = ('HIVTransTool',
                  'TreeingTools',
                  'HIVSeqDBManagement')

CELERY_MAX_CACHED_RESULTS = 100000
CELERY_TRACK_STARTED = True
CELERYD_MAX_TASKS_PER_CHILD = 10
CELERYD_POOL_RESTARTS = True
CELERY_SEND_EVENTS = True
CELERYD_HIJACK_ROOT_LOGGER = False
CELERYD_PREFETCH_MULTIPLIER = 1
CELERY_ACKS_LATE = True

CELERY_QUEUES = [Queue('HIVTransTool'),
                 Queue('celery'),
                 Queue('long-running'),
                 Queue('HIVSeqDBManagement')]

CELERY_ANNOTATIONS = {
   'HIVSeqDBManagement.query_LANL': {'rate_limit': '10/m'}
}