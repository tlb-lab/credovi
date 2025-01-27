#!/usr/bin/python

from credovi import app , run
from credovi.run import base, credo, preparepdb, mmcif, credosa

def main():
    """
    Dispatch the commands from here.
    """
    # run the application
    app.run()

    # get the controller that was specified on the command line
    controller = app.controller

    # this script was called without any controller argument
    if not controller: pass

    # interact with the current configuration
    elif controller.Meta.label == 'base': base.do(controller)

    # work with the CREDO schema
    elif controller.Meta.label == 'db':
        from credovi.run import db
        db.do(controller)

    elif controller.Meta.label == 'credo':

        # generate interaction data
        if controller.command == 'contacts': credo.do(controller)
        elif controller.command == 'preparepdb': preparepdb.do(controller)

    elif controller.Meta.label == 'mmcif': mmcif.do(controller)

    # run stand-alone version
    elif controller.Meta.label == 'credosa': credosa.do(controller)

    # run stand-alone version
    elif controller.Meta.label == 'ligand':

        if controller.command == 'fuzcav':
            from credovi.run import fuzcav
            fuzcav.do(controller)

        elif controller.command == 'molstrings':
            from credovi.run import molstrings
            molstrings.do(controller)

        elif controller.command == 'surfareas':
            from credovi.run import surfareas
            surfareas.do(controller)

    # close the application
    app.close()

if __name__ == '__main__':
    main()