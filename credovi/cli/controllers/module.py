from cement2.core import controller

class ModuleController(controller.CementBaseController):
    '''
    Each controller accepts its own commands, options and arguments.
    '''
    class Meta:
        interface = controller.IController
        label = 'module'
        description = """Run additional CREDO modules."""
        stacked_on = None
        config_defaults = dict(progressbar=False)
        arguments = [
            (['-P', '--progressbar'],
                dict(action='store_true', help='use a progressbar')),
            (['--dry-run'],
                dict(action='store_true', help='generate data but do not insert')),
            ]

    @controller.expose(hide=True, aliases=['run'])
    def default(self):
        '''
        Default method that is called if no commands were specified.
        '''
        print self.usage_text
        print
        print self.help_text

    @controller.expose(hide=False, help='create and insert new FuzCav fingerprints for ligands in CREDO')
    def fuzcav(self):
        '''
        '''
        pass