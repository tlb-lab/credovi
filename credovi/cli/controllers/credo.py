from cement.core import controller
from textwrap import fill

class CredoController(controller.CementBaseController):
    """
    Each controller accepts its own commands, options and arguments.
    """
    class Meta:
        interface = controller.IController
        label = 'credo'
        description = fill("Calculate the core CREDO interaction data for biological "
                           "complexes. Prepare PDB files for use with the CREDO database.", 80)
        stacked_on = None
        config_defaults = dict(testset='small', progressbar=False)

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
                     help='Only process structures that have not been processed yet')),
                     
            (['-U', '--update'],
                dict(nargs='?', metavar='DAYS', type=int, default=0, const=7,
                     help='Only process new structures or those whose output '
                          'files are older than DAYS (default: 7)')),

            (['-I', '--incremental'],
                dict(action='store_true',
                     help='Skip structures that have been processed before')),

            (['-T', '--testset'],
                dict(action='store', metavar='SIZE',
                     help='Only process PDB entries given in the test set')),

            (['--chains'],
                dict(action='store_true',
                     help="Save all chains in the prepared structure as separate files.")),

            (['--ligands'],
                dict(action='store_true',
                     help="Save all ligands in the prepared structure as separate files.")),

            (['--oeb'],
                dict(action='store_true',
                     help="Save the prepared structure in OpenEye binary format.")),

            (['--pdb'],
                dict(action='store_true',
                     help="Save the prepared structure in PDB format.")),

            (['--progressbar'],
                dict(action='store_true',
                     help="Display a progress bar and supress INFO messages.")),

            (['--clean'],
                dict(action='store_true',
                     help="Empty directory containing quaternary assemblies of a PDB entry first.")),

            (['--dry-run'],
                dict(action='store_true',
                     help="Process data but do not write files.")),
            ]

    @controller.expose(hide=True, aliases=['run'])
    def default(self):
        """
        Default method that is called if no commands were specified.
        """
        pass

    @controller.expose(hide=False, help="Default mode to generate the data files "
                       "for the CREDO database - requires access to local filesystems "
                       "and databases.")
    def contacts(self):
        """
        """
        pass

    @controller.expose(hide=False, help="Process PDB files from the PDB mirror. "
                       "Required for database generation.")
    def preparepdb(self):
        """
        python credovi.py mmcif currentpdbs -L 32 | parallel --dry-run --halt 2 -n 8 python credovi.py credo preparepdb -Q -S {1},{2},{3},{4},{5},{6},{7},{8}
        """
        pass

class CredoStandAloneController(controller.CementBaseController):
    """
    This controller is used to generate contacts on local files without any extra
    dependencies such as databases.
    """
    class Meta:
        interface = controller.IController
        label = 'credosa'
        description = """Calculate the core CREDO interaction data for biological complexes. Prepare PDB files for use with the CREDO database. Stand-alone version."""
        stacked_on = None
        config_defaults = dict(fmt='tsv')

        arguments = [
            (['-S', '--structures'],
                dict(action='store', metavar='PDB CODES',
                     help='PDB structures that should be processed')),
            (['-P', '--pymol'],
                dict(action='store_true',
                     help='Generate a PyMOL script to visualise contacts.')),
            (['-C', '--contacts'],
                dict(action='store_true',
                     help='Write contacts to a file.')),
            (['-F', '--fmt'],
                dict(action='store', metavar='FORMAT', default='tsv',
                     help="File format in which contacts will be saved in, either 'csv', 'tsv' or 'xls'. Defaults to 'tsv'.")),
            (['-O', '--output-dir'],
                dict(action='store',
                     help="Consider intramolecular interactions as well if possible.")),
            (['--gz'],
                dict(action='store_true',
                     help="Treat input files as being gzipped.")),
            (['--intramolecular'],
                dict(action='store_true',
                     help="Consider intramolecular interactions as well if possible. Otherwise only interactions between disconnected components will be recorded.")),
            (['--ri'],
                dict(action='store_true',
                     help="Show ring interactions.")),
            (['--undefined'],
                dict(action='store_true',
                     help="Record undefined contacts as well.")),
            (['--nolabels'],
                dict(action='store_true',
                     help="Hide distance labels.")),
            (['--progressbar'],
                dict(action='store_true',
                     help="Display a progress bar and supress INFO messages.")),
            (['--bindingsites'],
                dict(action='store_true',
                     help="Only consider protein-ligand binding sites for contact identification.")),
            (['--water'],
                dict(action='store_true',
                     help="Include water-mediated hydrogen bonds."))
            ]

    @controller.expose(hide=True, aliases=['run'])
    def default(self):
        """
        Default method that is called if no commands were specified.
        """
        pass