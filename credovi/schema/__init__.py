import json

from sqlalchemy import create_engine, DDL, MetaData
from sqlalchemy.event import listen

from credovi import app

__all__ = ['tables']

# engine depends on the configuration file that is parsed during app.setup()
engine = create_engine(app.config.get('database','url'),
                       echo=app.config.get('database','echo'))


schema = app.config.get('database','schema')
metadata = MetaData(bind=engine, schema=schema)

from .tables import *

listen(metadata, "after_create",
       DDL("GRANT SELECT ON ALL TABLES IN SCHEMA {} TO credouser".format(schema)))
