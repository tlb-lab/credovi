import os
import json
from optparse import OptionParser

from pychem import *
from openeye.oeomega import OEOmega
from progressbar import ProgressBar, Percentage, Bar, SimpleProgress

PDBCHEM_SDF_DIR  = '/tlbnas/mirror/pdbe/pdbechem/files/sdf_nh'
PDBCHEM_CONF_DIR = '/tlbnas/store/pdbe/pdbechem/conformers'

def parse_options():
    '''
    '''
    # PARSE COMMAND LINE
    # USE WITH PARALLEL: find /tlbnas/mirror/pdbe/pdbechem/files/sdf_nh -name *.sdf | awk -F"/" '{print $NF}' | parallel --eta --halt 2 -n 6 python generate-conformers.py -S -F \'{}\'
    usage  = "%prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-D", "--debug",
                      action  = "store_true",
                      dest    = "debug",
                      default = False,
                      help    = 'Set logging level to debug and print more verbose output.')

    parser.add_option("-S", "--silent",
                      action  = "store_true",
                      dest    = "silent",
                      default = True,
                      help    = 'Suppress Omega output.')

    parser.add_option("-C", "--chemcomps",
                      dest    = "chemcomps",
                      default = None,
                      help    = "PDB entries to include.")

    parser.add_option("-F", "--files",
                      dest    = "files",
                      default = None,
                      help    = "")

    parser.add_option("-I", "--incremental",
                      action  = "store_true",
                      dest    = "incremental",
                      default = False,
                      help    = "Only process chemical components that have not been processed before.")

    parser.add_option("-P", "--progressbar",
                      action  = "store_true",
                      dest    = "progressbar",
                      default = False,
                      help    = "Show a progress bar.",)

    parser.add_option("--surfaces",
                      action  = "store_true",
                      dest    = "surfaces",
                      default = False,
                      help    = '')

    # GET COMMAND LINE OPTIONS
    (options, args) = parser.parse_args()

    return options

def get_sdfs(options):
    '''
    '''
    sdfs = []

    # RETURN ALL SDFS IN DIRECTORY
    if options.files:
        for f in options.files.split(' '):
            sdf_path = os.path.join(PDBCHEM_SDF_DIR, f)
            if os.path.exists(sdf_path): sdfs.append(f)

    # INCLUDE ONLY THE SPECIFIED CHEMICAL COMPONENTS
    elif options.chemcomps:
        het_ids = set(het_id.strip().upper() for het_id in options.chemcomps.split(' '))

        for het_id in het_ids:
            sdf = '%s.sdf' % het_id
            sdf_path = os.path.join(PDBCHEM_SDF_DIR, sdf)
            if os.path.exists(sdf_path): sdfs.append(sdf)

    else:
        sdfs = os.listdir(PDBCHEM_SDF_DIR)

    return sorted(sdfs)

def main():
    '''
    '''
    options = parse_options()

    if options.silent: OEThrow.SetLevel(OEErrorLevel_Error)

    # DEFINE OMEGA PARAMETERS
    omega = OEOmega()
    omega.SetMaxConfs(20)
    omega.SetMaxRotors(30)
    omega.SetWarts(True)
    omega.SetSDEnergy(True)

    ifs = oemolistream()
    ifs.SetFormat(OEFormat_SDF)

    # WRITER FOR OUTPUT MOLECULES
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_OEB)

    sdfs = get_sdfs(options)

    # INITIALIZE PROGRESSBAR
    if options.progressbar:
        bar = ProgressBar(widgets=['Chemical components: ', SimpleProgress(), ' ', Percentage(), Bar()], maxval=len(sdfs)).start()

    for counter,sdf in enumerate(sdfs,1):
        if options.progressbar: bar.update(counter)
        else: print sdf

        het_id, ending = os.path.splitext(sdf)

        # IGNORE NON-SDF FILES
        if ending != '.sdf': continue

        sdf_path = os.path.join(PDBCHEM_SDF_DIR, sdf)
        conf_path = os.path.join(PDBCHEM_CONF_DIR, '%s.oeb' % het_id)

        # DO NOT REGENERATE CONFORMERS IF INCREMENTAL OPTION WAS CHOSEN
        if options.incremental and os.path.lexists(conf_path):
            continue

        # READ-IN SINGLE-CONFORMER SDF FILE
        ifs.open(str(sdf_path))

        # IGNORE DODGY MOLECULES
        try: chemcomp = ifs.GetOEMols().next()
        except StopIteration: continue

        chemcomp.SetTitle(str(het_id))

        # RESET CHARGES
        reset_charges(chemcomp)

        # CREATE A NEW MOLECULE FOR OMEGA BECAUSE IF OMEGA FAILS THE MOLECULE WILL
        # BE CLEARED
        chemcompconfs = OEMol(chemcomp)

        # GENERATE CONFORMERS WITH SPECIFIED PARAMETERS
        # FALL BACK TO SINGLE CONFORMER IF OMEGA FAILS
        if not omega(chemcompconfs): chemcompconfs = chemcomp

        # IGNORE HYDROGENS FOR SURFACE AREA CALCULATIONS
        OESuppressHydrogens(chemcompconfs)

        if options.surfaces:

            # CALCULATE SURFACE AREAS FOR EACH CONFORMER
            for conformer in chemcompconfs.GetConfs():
                asa,aasa,pasa = get_solvent_accessible_surface_areas(conformer)

                conformer.SetFloatData('asa',asa)
                conformer.SetFloatData('aasa',aasa)
                conformer.SetFloatData('pasa',pasa)

        # WRITE CONFORMERS
        ofs.open(str(conf_path))
        OEWriteMolecule(ofs, chemcompconfs)

    if options.progressbar: bar.finish()

main()
