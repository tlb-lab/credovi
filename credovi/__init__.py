import os

from sqlalchemy import create_engine
from cement2.core import backend, foundation, handler

from cli import JSONConfigParserHandler
from cli import controllers

# create an application
defaults = backend.defaults('credovi')

defaults['base']['config_handler'] = 'JSONConfigParser'
defaults['base']['config_files'] = [os.path.join(__path__[0], './config', 'config.json'),
                                    os.path.join(__path__[0], './config', 'credo.json'),
                                    os.path.join(__path__[0], './config', 'pymol.json')]

# defined app first and handlers afterwards
app = foundation.lay_cement('credovi', defaults=defaults)

# hooks must be registered after foundation.lay_cement()
from cli.hooks import *

handler.register(JSONConfigParserHandler)
handler.register(controllers.BaseController)
handler.register(controllers.DatabaseController)
handler.register(controllers.CredoController)
handler.register(controllers.CredoStandAloneController)
handler.register(controllers.MMCIFController)

# setup the application
app.setup()

# engine depends on the configuration file that is parsed during app.setup()
engine = create_engine(app.config.get('database','url'),
                       echo=app.config.get('database','echo'))