"""
This module contains functions to fetch structural data from a database, most
notably the mmcif database that contains the complete mmCIF data.
"""
import time
from itertools import groupby
from operator import itemgetter
from sqlalchemy.sql import text

from openeye.oechem import OEFloatArray
from credovi import app
from credovi.schema import engine

def get_credo_pdbs(generator=True):
    """
    Returns a list of all the PDB entries that are currently stored in CREDO.
    Important for doing incremental updates.
    """
    statement = text("SELECT pdb FROM credo.structures ORDER BY 1")
    result = engine.execute(statement).fetchall()

    if generator: return (row.pdb for row in result)

    return result

def get_pdb_polymer_info(pdb):
    """
    Returns a complete polymer sequence residue mapping including the entity type
    bit mask.
    """
    statement = text("""
                    SELECT      pdb_strand_id as pdb_chain_id, pdb_mon_id as res_name,
                                pdb_seq_num::integer AS res_num,
                                CASE
                                    WHEN pdb_ins_code = '' THEN ' '
                                    ELSE pdb_ins_code
                                END AS ins_code,
                                CASE
                                    -- PROTEIN
                                    WHEN ep.type = 'polypeptide(L)' OR ep.type = 'polypeptide(D)' THEN 32
                                    -- DNA
                                    WHEN ep.type = 'polydeoxyribonucleotide' THEN 16
                                    -- RNA
                                    WHEN ep.type = 'polyribonucleotide' THEN 8
                                    -- DNA/RNA HYBRID
                                    WHEN ep.type = 'polydeoxyribonucleotide/polyribonucleotide hybrid' THEN 24
                                    -- POLYSACCHARIDE
                                    WHEN ep.type = 'polysaccharide(D)' THEN 4
                                    ELSE 0
                                END AS entity_type_bm
                    FROM        mmcif.pdbx_poly_seq_scheme p
                    JOIN        mmcif.entity_poly ep ON p.structure_id = ep.structure_id AND p.entity_id = ep.entity_id
                    WHERE       p.pdb_mon_id != '' AND p.Structure_ID = :pdb
                    ORDER BY    1, 3, 4
                     """)

    result = engine.execute(statement, pdb=pdb.upper()).fetchall()

    # create a mapping between the pdb identifier of the residue and its entity type
    residues = (((row.pdb_chain_id, row.res_name, row.res_num, row.ins_code), row.entity_type_bm) for row in result)
    residues = dict(residues)

    app.log.debug("structure contains {0} polymer residues according to mmCIF.".format(len(residues)))

    return residues

def get_pdb_ligand_info(pdb):
    """
    """
    SQL = text("""
                SELECT      pdb_chain_id, het_id, res_num, ins_code, 2 as entity_type_bm
                FROM        pdb.ligands
                WHERE       pdb = :pdb AND ins_code = ' ' AND is_observed = true
                UNION
                SELECT      pdb_chain_id, NULL, NULL, ' ', 34 as entity_type_bm
                FROM        pdb.peptide_ligands
                WHERE       pdb = :pdb
                ORDER BY    1, 3, 4
               """)

    result = engine.execute(SQL, pdb=pdb.upper()).fetchall()

    # create a mapping between the pdb identifier of the residue and its entity type
    ligands = (((row.pdb_chain_id, row.het_id, row.res_num, row.ins_code), row.entity_type_bm) for row in result)
    ligands = dict(ligands)

    app.log.debug("structure contains {0} ligands according to mmCIF:"
                  .format(len(ligands)))

    for tup in ligands.keys():
        app.log.debug("    {0} {1} {2}".format(*tup))

    return ligands

def get_pdb_sstruct_info(pdb):
    """
    """
    statement = text("""
                    SELECT  pdb_chain_id, pdb_res_name, pdb_res_num, pdb_ins_code, sstruct_serial
                    FROM    pdb.res_map
                    WHERE   pdb = :pdb
                    """)

    result = engine.execute(statement, pdb=pdb.upper()).fetchall()
    sstruct_info = {}

    # iterate through secondary structure mapping for this pdb sequence
    for row in result:

        # list is not hashable
        row = tuple(row.values())

        # create a key/value pair in the form pdb id => sstruct serial
        sstruct_info[row[:-1]] = row[-1]

    app.log.debug("{0} protein fragments were identified through SIFTS."
                  .format(len(sstruct_info)))

    return sstruct_info

def get_pisa_data(pdb):
    """
    """
    statement = text("""
                SELECT      s.serial as set_serial, a.serial as assembly_serial, a.mmsize as assembly_size, pdb_select, rotation, translation, is_at_identity
                FROM        pisa.entries e
                JOIN        pisa.sets s ON s.entry_id = e.entry_id
                JOIN        pisa.assemblies a ON a.set_id = s.set_id
                JOIN        pisa.molecules m ON m.assembly_id = a.assembly_id
                WHERE       -- CHECK IF PISA PREDICTION EXISTS
                            e.status = 'Ok'
                            -- FIRST SET IS ALWAYS THE TOP-RATED
                            AND s.serial = 1
                            -- ONLY STABLE ASSEMBLIES
                            AND a.is_stable = true
                            -- ONLY ASSEMBLIES CONTAINING BIOMOLECULES
                            AND a.has_polymer = true
                            -- FOR PDB CODE
                            AND e.pdb = :pdb
                            -- ONLY COMPLETE CHAINS
                            AND LENGTH(pdb_select) = 1
                ORDER BY    s.serial, a.serial;
                """)

    # fetch data from pisa database
    result = engine.execute(statement, pdb=pdb.upper()).fetchall()

    pisa = {}
    ASSEMBLY_TOO_LARGE = False

    # transform result into dictionary data structure
    for set_serial, assembly_iter in groupby(result, key=itemgetter(0)):
        pisa[set_serial] = {}

        # iterate through assemblies
        for (assembly_serial, assembly_size), chain_iter in groupby(assembly_iter, key=itemgetter(1,2)):
            pisa[set_serial].update({assembly_serial:{}})

            # warn if an assembly will be ignored
            # only 194 pdb structures have assemblies with more than 26 chains
            if assembly_size > 26: ASSEMBLY_TOO_LARGE = True

            # iterate through chains
            for pdb_chain_id, operation_iter in groupby(chain_iter, key=itemgetter(3)):
                pisa[set_serial][assembly_serial].update({str(pdb_chain_id):{}})

                for operation_serial, operation in enumerate(operation_iter,1):
                    rotation, translation, is_at_identity = operation[4:]

                    details = {'rotation': OEFloatArray(rotation),
                               'translation': OEFloatArray(translation),
                               'is_at_identity': is_at_identity}

                    pisa[set_serial][assembly_serial][pdb_chain_id][operation_serial] = details

    # debug assembly information
    app.log.debug("PISA predicts {0} assembly/assemblies."
                  .format(max(pisa[1].keys())))

    # do not return pisa data for large assemblies / asu will be used instead
    if ASSEMBLY_TOO_LARGE:
        app.log.warn("At least one of the predicted assemblies is too large - "
                     "the asymmetric unit will be used instead.")
        pisa = {}

    return pisa

def get_pisa_num_assembly_sets(pdb):
    """
    """
    statement = text("""
                          SELECT num_assembly_sets, s.serial, a.is_stable, a.has_polymer
                            FROM pisa.entries e
                       LEFT JOIN pisa.sets s ON s.entry_id = e.entry_id
                       LEFT JOIN pisa.assemblies a ON a.set_id = s.set_id
                           WHERE pdb = :pdb AND status = 'Ok'
                        ORDER BY 2,3,4;
                     """)

    result = engine.execute(statement, pdb=pdb.upper()).fetchall()

    for row in result:
        num_assembly_sets, set_serial, is_stable, has_polymer = row.values()

        # PISA predicts monomer, ASU has to be split
        if num_assembly_sets == 0: return 0

        # return the number of assemblies if stable
        if set_serial == 1 and is_stable and has_polymer: return num_assembly_sets

    # no PISA entry, proceed with ASU
    return -1
