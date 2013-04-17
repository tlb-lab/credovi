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
            (['-A', '--all'],
                dict(action='store_true',
                     help='Process all structures that can be found in the appropriate directory for the given task.')),

            (['-Q', '--quat'],
                dict(action='store_true',
                     help='Only consider quaternary assemblies')),

            (['-S', '--structures'],
                dict(action='store', metavar='PDB CODES',
                     help='PDB structures that should be processed')),

            (['-O', '--offset'],
                dict(action='store', metavar='PDB CODE',
                     help='PDB code to use as an offset for processing')),

            (['-N', '--new'],
                dict(action='store_true',
                     help='Only structures that have not been processed yet')),

            (['-I', '--incremental'],
                dict(action='store_true',
                     help='only process new data')),

            (['-T', '--testset'],
                dict(action='store', metavar='SIZE',
                     help='Only process PDB entries given in the test set')),

            (['-P', '--progressbar'],
                dict(action='store_true', help='use a progressbar')),
            (['--dry-run'],
                dict(action='store_true', help='generate data but do not insert')),
            (['--usr'],
                dict(action='store_true', help='also generate and insert USR moments for ligand structures')),
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

    @controller.expose(hide=False, help="calculate buried surface areas")
    def surfareas(self):
        """
        """
        pass