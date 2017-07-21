from cement.core import controller

class BaseController(controller.CementBaseController):
    """
    Base controller for Credovi. This controller will be called if no other controller
    was specified; other controllers can be stacked on this one if necessary.
    """
    class Meta:
        """
        """
        interface = controller.IController
        label = 'base' # do not name 'base' otherwise it will clash with config
        description = "Credovi: command line interface to work with the CREDO database."

        # default values for command line options
        config_defaults = dict()

        # definition of command line options that are supported
        arguments = [
            (['-S', '--section'],
                dict(action='store',
                     help='Print all keys of the given section.')),
            (['-K', '--keys'],
                dict(action='store',
                     help='Print all keys with the given name of only of a specific section.')),
            (['--single'],
                dict(action='store_true',
                     help='Only print one value per line.')),
            (['--sort'],
                dict(action='store_true',
                     help='Sort return values if possible.')),
            (['--list-not-in-db'],
                dict(action='store_true',
                     help='List the PDB entries that are not stored in CREDO yet.')),
            ]

    @controller.expose(hide=True, aliases=['run'])
    def default(self):
        """
        """
        print "FOO"

    @controller.expose(help="Print sections and keys of the current configuration.")
    def conf(self):
        """
        Do not name config, otherwise it will clash with other objects attached
        to the app.
        """
        pass

    @controller.expose(help="Work with the CREDO data directory.")
    def data(self):
        """
        """
        pass