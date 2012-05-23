from cement2.core import controller

class MMCIFController(controller.CementBaseController):
    '''
    Each controller accepts its own commands, options and arguments.
    '''
    class Meta:
        interface = controller.IController
        label = 'mmcif'
        description = "Interact with the mmcif database."
        stacked_on = None
        config_defaults = dict(limit=5000)

        arguments = [
            (['-O', '--offset'], dict(action='store', metavar='PDB CODE', help='PDB code to use as an offset for processing.')),
            (['-L', '--limit'], dict(action='store', metavar='NUMBER', help="Maximum number of PDB codes to return.")),
            ]

    @controller.expose(hide=True, aliases=['run'])
    def default(self):
        '''
        Default method that is called if no commands were specified.
        '''
        pass

    @controller.expose(hide=False, help='Get the (sorted) PDB codes of all the structures that are currently stored in the mmcif database. Normally used in conjunction with parallel.')
    def currentpdbs(self):
        '''
        '''
        pass