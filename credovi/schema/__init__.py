import json

from sqlalchemy import MetaData

from credovi import app, engine

__all__ = ['tables']

metadata = MetaData(bind=engine)
schema = app.config.get('database','schema')

from .tables import *