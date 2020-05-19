from __future__ import absolute_import, unicode_literals
from celery import Celery

BROKER = 'amqp://myuser:mypassword@localhost:5672/myvhost'
BACKEND = 'redis://localhost:6379/0'

app = Celery('celery1', backend=BACKEND, broker=BROKER,include=['celsol.tasks'])
app.conf.update(
    broker_transport_options = {'visibility_timeout': 3600} # 1 h
)

if __name__ == '__main__':
    app.start()

