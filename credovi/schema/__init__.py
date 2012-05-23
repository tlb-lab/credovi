import json

from sqlalchemy import create_engine, MetaData

from credovi import app

__all__ = ['tables']

# engine depends on the configuration file that is parsed during app.setup()
engine = create_engine(app.config.get('database','url'),
                       echo=app.config.get('database','echo'))

metadata = MetaData(bind=engine)
schema = app.config.get('database','schema')

from .tables import *