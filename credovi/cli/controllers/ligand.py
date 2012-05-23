from cement2.core import controller

class LigandController(controller.CementBaseController):
    """
    """
    class Meta:
        interface = controller.IController
        label = 'ligand'
        description = """Run additional CREDO modules."""
        stacked_on = None
        config_defaults = dict(progressbar=False)
        arguments = [
            (['-P', '--progressbar'],
                dict(action='store_true', help='use a progressbar')),
            (['--dry-run'],
                dict(action='store_true', help='generate data but do not insert')),
            (['--usr'],
                dict(action='store_true', help='also generate and insert USR moments for ligand structures')),
            (['--incremental'],
                dict(action='store_true', help='only process new data')),
            ]

    @controller.expose(hide=True, aliases=['run'])
    def default(self):
        """
        Default method that is called if no commands were specified.
        """
        print self.usage_text
        print
        print self.help_text

    @controller.expose(hide=False, help="create and insert new FuzCav fingerprints for ligands in CREDO")
    def fuzcav(self):
        """
        """
        pass

    @controller.expose(hide=False, help="insert ligand structures in different file formats")
    def molstrings(self):
        """
        """
        pass