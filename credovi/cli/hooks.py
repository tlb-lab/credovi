import sys

from cement2.core import hook
from openeye.oechem import OEChemIsLicensed

@hook.register(name='cement_pre_run_hook')
def openeye_is_licensed_hook(app):
    """
    Check if a valid OpenEye license is present before running the application.
    """
    if hasattr(app.controller, 'Meta'):
        if app.controller.Meta.label == 'credo' and not OEChemIsLicensed():
            app.log.fatal("a valid OpenEye license is required for running credovi "
                          "with the 'credo' controller.")
            app.close()

@hook.register(name='cement_post_run_hook')
def check_cmd_args_hook(app):
    """
    Check the combination of command line arguments that have been provided. Some
    are mutually exclusive and testing them through a hook seems to be the best
    solution.
    """
    controller = app.controller
    
    # no controller, nothing to check
    if not controller: return
    
    # get the controller command
    cmd = controller.command
    
    # get the command line arguments and options
    args = controller.pargs    
    
    # debugging and displaying a progress bar does not work well obviously
    if args.debug and args.progressbar:
        app.log.fatal("The --debug flag cannot be combined with --progressbar.")
        app.close()
    
    if controller.Meta.label == 'credo':
        if args.all:
            
            # --all flag cannot be combined with some of the other flags
            if any((args.structures, args.new, args.testset)):
                app.log.fatal("The --all flag cannot be combined with --structures, --new or --testset.")
                app.close()
        
        elif args.testset:
        
            # ditto for --testset
            if any((args.structures, args.new)):
                app.log.fatal("The --testset flag cannot be combined with --structures or --new.")
                app.close()
            
            elif args.testset not in ('small', 'large'):
                app.log.fatal("--testset must be either 'small' or 'large'.")
                app.close()
    
    elif controller.Meta.label == 'credosa':
        if not args.structures:
            app.log.fatal("a list of files must be provided with the --structures option.")
            app.close()
            sys.exit(1)