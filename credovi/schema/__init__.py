import json
from getpass import getpass

from sqlalchemy import create_engine, DDL, MetaData
from sqlalchemy.event import listen
from sqlalchemy.engine.url import make_url

from credovi import app

__all__ = ['tables']

# engine depends on the configuration file that is parsed during app.setup()
db_url = make_url(app.config.get('database','url'))
if not db_url.username:
    db_url.username = app.config.get('database','username')
if not db_url.password:
    db_url.password = getpass("Input password for user %s: " % db_url.username) if 'password' not in app.config.keys('database') \
                      else app.config.get('database','password')


engine = create_engine(db_url,
                       echo=app.config.get('database','echo'))


schema = app.config.get('database','schema')
metadata = MetaData(bind=engine, schema=schema)

from .tables import *

listen(metadata, "after_create",
       DDL("GRANT SELECT ON ALL TABLES IN SCHEMA {} TO credouser".format(schema)))
