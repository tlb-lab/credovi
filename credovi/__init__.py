import os
import glob

from cement2.core import backend, foundation, handler

from cli import JSONConfigParserHandler
from cli import controllers

CONFIG_DIR = os.path.join(__path__[0], './config')
CONFIG_FILES = glob.glob(os.path.join(CONFIG_DIR, '*.json'))

# defined app first and handlers afterwards
app = foundation.CementApp('credovi',
                           config_defaults=backend.defaults('credovi'),
                           config_handler='JSONConfigParser',
                           config_files=CONFIG_FILES)

# hooks must be registered afterwards
from cli.hooks import *

handler.register(JSONConfigParserHandler)
handler.register(controllers.BaseController)
handler.register(controllers.DatabaseController)
handler.register(controllers.CredoController)
handler.register(controllers.CredoStandAloneController)
handler.register(controllers.MMCIFController)
handler.register(controllers.LigandController)

handler.register(controllers.ModuleController)

# setup the application
app.setup()