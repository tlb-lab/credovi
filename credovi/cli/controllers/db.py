from textwrap import fill
from cement.core import controller

class DatabaseController(controller.CementBaseController):
    """
    Each controller accepts its own commands, options and arguments.
    """
    class Meta:
        interface = controller.IController
        label = 'db'
        description = fill("Manage the CREDO database schema as well as the data. "
                           "This controller can be used to CREATE/DROP/TRUNCATE tables "
                           "in CREDO but also to load new data files into the raw "
                           "tables for subsequent insertion into CREDO.", 80)
        stacked_on = None
        config_defaults = dict()
        arguments = [
            (['-t', '--tables'], dict(action='store', help='only the specified tables - if empty, ALL currently defined tables will be used')),
            (['-cf', '--checkfirst'], dict(action='store_true', help='check first (SQLAlchemy MetaData) before dropping/creating a table')),
            (['--echo'], dict(action='store_true', help='echo all issued SQL statements')),
            (['--raw'], dict(action='store_true', help="only include the 'raw' tables that are used for loading data")),
            (['--core'], dict(action='store_true', help="only include the 'core' tables")),
            (['--missing'], dict(action='store_true', help="only include missing tables (i.e. do not drop pre-existing) (creation only)")),
            (['--cascade'], dict(action='store_true', help='CASCADE all drop statements')),
            (['--sure'], dict(action='store_true', help='extra flag that needs to be specified for certain operations on the database')),
            (['--new'], dict(action='store_true', help='copy the data files of all PDB entries that are not in CREDO yet')),
            (['--no-indexes'], dict(action='store_true', help='')),
            (['--no-constraints'], dict(action='store_true', help='')),
            (['--add-depend'], dict(action='store_true', help='Also add dependent tables')),
            (['--parallel'], dict(action='store_true', help='use GNU parallel for certain tasks (cryst only!)')),
            ]

    @controller.expose(hide=True, aliases=['run'])
    def default(self):
        """
        Default method that is called if no commands were specified.
        """
        pass

    @controller.expose(hide=False, help='copy new data from the CREDO data directory into the raw tables in the database')
    def copy(self):
        """
        """
        pass

    @controller.expose(hide=False, help='CREATE all currently defined schema elements')
    def create(self):
        """
        """
        pass

    @controller.expose(hide=False, help='DROP all tables in the CREDO database schema')
    def drop(self):
        """
        """
        pass

    @controller.expose(hide=False, help='print the DDL of currently defined tables')
    def dump(self):
        """
        """
        pass

    @controller.expose(hide=False, help='TRUNCATE tables in the CREDO database schema')
    def truncate(self):
        self.log.info('Inside controller3.command3 function.')