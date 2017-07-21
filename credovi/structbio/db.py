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
                    FROM        {mmcif}.pdbx_poly_seq_scheme p
                    JOIN        {mmcif}.entity_poly ep ON p.structure_id = ep.structure_id AND p.entity_id = ep.entity_id
                    WHERE       p.pdb_mon_id != '' AND p.Structure_ID = :pdb
                    ORDER BY    1, 3, 4
                     """.format(mmcif=app.config.get('database','mmcif')))

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
                FROM        pdb_dev.ligands
                WHERE       pdb = :pdb AND (ins_code = '' OR ins_code = ' ') -- AND is_observed = true
                UNION
                SELECT      pdb_chain_id, NULL, NULL, ' ', 34 as entity_type_bm
                FROM        pdb_dev.peptide_ligands
                WHERE       pdb = :pdb
                ORDER BY    1, 3, 4
               """)

    result = engine.execute(SQL, pdb=pdb.upper()).fetchall()

    # create a mapping between the pdb identifier of the residue and its entity type
    ligands = ( ((row.pdb_chain_id, row.het_id, row.res_num, row.ins_code or ' '), row.entity_type_bm) for row in result)
    ligands = dict(ligands)

    app.log.debug("structure {0} contains {1} ligands according to mmCIF:"
                  .format(pdb, len(ligands)))

    for tup in ligands.keys():
        app.log.debug("    {0} {1} {2}".format(*tup))

    return ligands

def get_pdb_sstruct_info(pdb):
    """
    """
    statement = text("""
                    SELECT  pdb_chain_id, pdb_res_name, pdb_res_num, pdb_ins_code, sstruct_serial
                    FROM    pdb_dev.res_map
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

def get_biomt(pdb):
    """
    """
    statement = text("""
                        SELECT assembly_serial, assembly_size, is_monomeric,
                               pdb_chain_id, rotation, translation, is_at_identity
                          FROM pdb_dev.biomt b
                          JOIN pdb_dev.biomt_ops o USING(biomt_id)
                         WHERE pdb = :pdb
                      ORDER BY 1,3
                     """)

    # fetch data from pisa database
    result = engine.execute(statement, pdb=pdb.upper()).fetchall()

    biomt = {}
    is_monomeric = False

    # iterate through assemblies
    for (assembly_serial, assembly_size, is_monomeric), chain_iter in groupby(result, key=itemgetter(0,1,2)):
        biomt[assembly_serial] = {}

        # do not return pisa data for large assemblies / asu will be used instead
        if assembly_size > 26:
            app.log.warn("one of the predicted assemblies contains {} chains - "
                         "the asymmetric unit will be used instead."
                         .format(assembly_size))
            biomt = {}
            break

        # the complete ASU is monomeric and has to be split into individual chains
        if is_monomeric: break

        # iterate through chains
        for pdb_chain_id, operation_iter in groupby(chain_iter, key=itemgetter(3)):
            biomt[assembly_serial].update({str(pdb_chain_id):{}})

            for operation_serial, operation in enumerate(operation_iter,1):
                rotation, translation, is_at_identity = operation[4:]

                details = {'rotation': OEFloatArray(rotation),
                           'translation': OEFloatArray(translation),
                           'is_at_identity': is_at_identity}

                biomt[assembly_serial][pdb_chain_id][operation_serial] = details

    # debug assembly information
    try:
        app.log.debug("BIOMT contains {0} assembly/assemblies."
                      .format(max(biomt.keys())))

    # PDB entry does not have a REMARK 350
    except ValueError:
        app.log.debug("NO REMARK 350 found.")

    return biomt, is_monomeric
