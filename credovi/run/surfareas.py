"""
python credovi.py mmcif currentpdbs | parallel --eta --halt 2 -n 6 python credovi.py ligand surfareas -Q -S {1},{2},{3},{4},{5},{6}
"""
import os
import sys
import csv
from itertools import compress

import numpy as np
from openeye import oechem, oespicoli
from progressbar import ProgressBar, Percentage, Bar, SimpleProgress

from credovi import app
from credovi.schema import engine
from credovi.structbio.structure import get_assembly_sets

# these custom predicates are currently in the eyesopen module
from eyesopen.predicates import OEAtomHasIntData, OEAtomBinaryAndIntData

# register new dialect to write tab-delimited files
csv.register_dialect('tabs', delimiter='\t')

def insert():
    """
    """
    statement = """
                 DO $$
                    DECLARE
                        biomol_id INTEGER;
                    BEGIN
                        FOR biomol_id IN SELECT DISTINCT biomolecule_id FROM credo.ligands ORDER BY 1
                        LOOP
                            EXECUTE
                            '
                               INSERT INTO credo.binding_site_atom_surface_areas
                               SELECT l.ligand_id, a.atom_id, rw.asa_apo, rw.asa_bound, rw.asa_delta
                                 FROM credo.raw_binding_site_atom_surface_areas rw
                                 JOIN credo.structures s ON s.pdb = rw.pdb
                                 JOIN credo.biomolecules b
                                      ON b.structure_id = s.structure_id
                                      AND b.assembly_serial = rw.assembly_serial
                                 JOIN credo.ligands l
                                      ON l.biomolecule_id = b.biomolecule_id
                                      AND l.entity_serial = rw.entity_serial
                                 JOIN credo.atoms a
                                      ON a.biomolecule_id = b.biomolecule_id
                                      AND a.atom_serial = rw.atom_serial
                                WHERE a.biomolecule_id = $1
                             ORDER BY 1,2
                            ' USING biomol_id;
                            RAISE NOTICE 'inserted binding site atom surface areas for biomolecule %', biomol_id;
                        END LOOP;
                END$$;
                """

    connection.execute(statement)

def get_binding_site(structure, ligand):
    """
    """
    partlist = oechem.OEIntArray(structure.GetMaxAtomIdx())

    # get all contacts between ligand and assembly
    contacts = oechem.OEGetNearestNbrs(structure, ligand, 5.0)
    binding_site_atom_idxs = set(contact.GetBgn().GetIdx() for contact in contacts)

    # create the partition map
    for atom in structure.GetAtoms():
        if atom.GetIdx() in binding_site_atom_idxs:
            partlist[atom.GetIdx()] = 1

    entitypred = oechem.OEPartPredAtom(partlist)

    # select the binding site atoms
    entitypred.SelectPart(1)

    # create a new molecule for the entity
    binding_site = oechem.OEGraphMol()
    oechem.OESubsetMol(binding_site, structure, entitypred, False, False)

    return binding_site

def get_atom_surface_areas(molecule, surface):
    """
    """
    # create empty array to hold area values for all triangles
    areas = oechem.OEFloatArray(surface.GetNumTriangles())

    # fill array with calculated surface areas
    oespicoli.OECalculateTriangleAreas(surface, areas)

    # create empty array to hold atom surface areas
    atom_areas = oechem.OEFloatArray(molecule.GetMaxAtomIdx())

    # fill array with surface contributions for each atom
    for i in xrange(surface.GetNumTriangles()):

        # get the triangle elements
        v1 = surface.GetTrianglesElement(i * 3)
        v2 = surface.GetTrianglesElement(i * 3 + 1)
        v3 = surface.GetTrianglesElement(i * 3 + 2)

        # get the atom indices for each triangle element
        a1 = surface.GetAtomsElement(v1)
        a2 = surface.GetAtomsElement(v2)
        a3 = surface.GetAtomsElement(v3)

        atom_areas[a1] += areas[i] / 3.0
        atom_areas[a2] += areas[i] / 3.0
        atom_areas[a3] += areas[i] / 3.0

    return np.array(atom_areas)

def write_atoms(writer, molecule, atom_map, pdb, assembly_serial, entity_serial,
                atom_areas_apo, atom_areas_bound):
    """
    """
    for atom in compress(molecule.GetAtoms(), atom_map):
        asa_apo = atom_areas_apo[atom.GetIdx()]
        asa_bound = atom_areas_bound[atom.GetIdx()]
        asa_delta = asa_apo - asa_bound

        residue = oechem.OEAtomGetResidue(atom)

        row = ['\N' for i in range(7)]

        row[0] = pdb
        row[1] = assembly_serial
        row[2] = entity_serial
        row[3] = residue.GetSerialNumber()
        row[4]= asa_apo
        row[5]= asa_bound
        row[6]= asa_delta

        writer.writerow(row)

def do(controller):
    """
    """
    # get the controller command
    cmd = controller.command

    # get the command line arguments and options
    args = controller.pargs

    # predicate to remove non-polymer atoms from structure
    nonpolymers = oechem.OEOrAtom(OEAtomHasIntData(('entity_type_bm', 0)),
                                  OEAtomBinaryAndIntData(('entity_type_bm', 3)))

    assemblysets = get_assembly_sets(args)

    # directory containing all the biological assemblies in OEB format
    OEB_ASSEMBLIES_DIR = app.config.get('directories','quat_oeb')

    # directory where surface areas will be written
    CREDO_DATA_DIR = app.config.get('directories','credo_data')

    ifs = oechem.oemolistream()
    ifs.SetFormat(oechem.OEFormat_OEB)

    # initialize progressbar
    if args.progressbar:
        bar = ProgressBar(widgets=['PDB entries: ', SimpleProgress(), ' ',
                                   Percentage(), Bar()],
                          maxval=len(assemblysets)).start()

    # iterate through assembly sets
    for counter, (pdb, assemblyset) in enumerate(assemblysets, 1):
        if args.progressbar: bar.update(counter)

        # create a data directory for this structure to which all data will be written
        struct_data_dir = os.path.join(CREDO_DATA_DIR,
                                       pdb[1:3].lower(), pdb.lower())

        # make necessary directories recursively if they do not exist yet
        if not os.path.exists(struct_data_dir):
            os.makedirs(struct_data_dir)

        # path to the file where the atom surface areas of all atoms will be written
        surface_areas_path = os.path.join(struct_data_dir,
                                          'binding_site_atom_surface_areas.credo')

        # do not recalculate atom surface area contributions if incremental
        if args.incremental and os.path.exists(surface_areas_path):
            continue

        # output file stream and CSV writer
        atomfs = open(surface_areas_path, 'w')
        atomwriter = csv.writer(atomfs, dialect='tabs')

        # deal with each found assembly separately
        # some pdb entries consist of more than one
        for assembly in assemblyset:
            if args.quat:
                path = os.path.join(OEB_ASSEMBLIES_DIR, pdb[1:3].lower(),
                                    pdb.lower(), assembly)

            else:
                app.log.error("the calculation of buried ligand surface areas "
                              "is only supported for quaternary structures.")
                sys.exit(1)

            if not os.path.isfile(path):
                app.log.warn("cannot calculate buried surface areas: "
                             "file {} does not exist!".format(path))

            # get the quaternary structure
            ifs.open(str(path))
            assembly = ifs.GetOEGraphMols().next()

            if not assembly:
                app.log.warn("cannot calculate buried surface areas: "
                             "file {} does not contain a valid molecule!"
                             .format(path))
                continue

            if not assembly.GetListData('ligands'):
                continue

            # identifier of the assembly
            assembly_serial = assembly.GetIntData('assembly_serial')

            # remove all non-polymers from assembly
            for atom in assembly.GetAtoms(nonpolymers):
                assembly.DeleteAtom(atom)

            # ignore bizarre assemblies
            if not assembly.NumAtoms():
                app.log.warn("cannot calculate buried surface areas: "
                             "file {} contains assembly with no atoms!"
                             .format(path))
                continue

            # keep only the location state with the largest average occupancy
            assembly_hi_occ = oechem.OEGraphMol()
            altlocfactory = oechem.OEAltLocationFactory(assembly)
            altlocfactory.MakeCurrentAltMol(assembly_hi_occ)

            # get the ligands
            ligands = assembly_hi_occ.GetListData('ligands')

            # iterate through all ligands of the biomolecule and calculate the buried
            # surface area atom contributions for all involved atoms
            for ligand in ligands:

                # ignore small ligands
                if oechem.OECount(ligand, oechem.OEIsHeavy()) < 7: continue

                entity_serial = ligand.GetIntData('entity_serial')

                # keep only the location state with the largest average occupancy
                altlig = oechem.OEGraphMol()
                altlocfactory = oechem.OEAltLocationFactory(ligand)
                altlocfactory.MakeCurrentAltMol(altlig)

                cmplx_srf = oespicoli.OESurface()
                ligand_srf = oespicoli.OESurface()

                # make solvent-accessible surface of ligand
                oespicoli.OEMakeAccessibleSurface(ligand_srf, altlig, 0.5, 1.4)

                # get the atom contributions of the assembly surface
                ligand_atom_areas = get_atom_surface_areas(altlig, ligand_srf)

                # extract the binding site of the assembly to speed up surface
                # area calculation
                binding_site = get_binding_site(assembly_hi_occ, altlig)

                # make solvent-accessible surface of binding site
                binding_site_srf = oespicoli.OESurface()
                oespicoli.OEMakeAccessibleSurface(binding_site_srf,
                                                  binding_site, 0.5, 1.4)

                # get the atom contributions of the assembly surface
                binding_site_atom_areas = get_atom_surface_areas(binding_site,
                                                                 binding_site_srf)

                # create complex
                cmplx = oechem.OEGraphMol()
                oechem.OEAddMols(cmplx, binding_site)
                oechem.OEAddMols(cmplx, altlig)

                # make solvent-accessible surface of the complex
                oespicoli.OEMakeAccessibleSurface(cmplx_srf, cmplx, 0.5, 1.4)

                # surface area atom contributions of the whole complex
                cmplx_atom_areas = get_atom_surface_areas(cmplx, cmplx_srf)

                ## extract the atom surface areas in the bound state through slices
                binding_site_atom_areas_bound = cmplx_atom_areas[:binding_site.NumAtoms()]
                ligand_atom_areas_bound = cmplx_atom_areas[binding_site.NumAtoms():]

                # difference between apo and bound state per polymer atom
                binding_site_delta = binding_site_atom_areas - binding_site_atom_areas_bound
                ligand_delta = ligand_atom_areas - ligand_atom_areas_bound

                # boolean map indicating for which atom the surface area has changed
                binding_site_atom_map = binding_site_delta != 0
                ligand_atom_map = ligand_delta != 0

                if args.dry_run: continue

                # only record the atoms where the solvent-accessible surface
                # area has actually changed
                write_atoms(atomwriter, binding_site, binding_site_atom_map, pdb,
                            assembly_serial, entity_serial, binding_site_atom_areas,
                            binding_site_atom_areas_bound)

                # only record the atoms where the solvent-accessible surface area
                # has actually changed
                write_atoms(atomwriter, altlig, ligand_atom_map, pdb, assembly_serial,
                            entity_serial, ligand_atom_areas, ligand_atom_areas_bound)

                app.log.debug("wrote buried surface areas for all ligands in "
                              "biomolecule {} to {}.".format(pdb, surface_areas_path))

            atomfs.flush()
        atomfs.close()

    if args.progressbar:
        bar.finish()