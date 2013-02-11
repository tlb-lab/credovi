"""
"""

from sqlalchemy import Boolean, Column, DefaultClause, Float, Index, Integer, Table
from sqlalchemy.schema import PrimaryKeyConstraint

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, comment_on_table_elements

interfaces = Table('interfaces', metadata,
                   Column('interface_id', Integer, nullable=False),
                   Column('biomolecule_id', Integer, nullable=False),
                   Column('chain_bgn_id', Integer, nullable=False),
                   Column('chain_end_id', Integer, nullable=False),
                   Column('path', PTree),
                   Column('num_res_bgn', Integer, nullable=False),
                   Column('num_res_end', Integer, nullable=False),
                   Column('is_quaternary', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                   Column('has_incomplete_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                   Column('has_non_std_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                   Column('has_mod_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                   Column('has_mut_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                   schema=schema)

PrimaryKeyConstraint(interfaces.c.interface_id, deferrable=True, initially='deferred')
Index('idx_interfaces_chain_bgn_id', interfaces.c.chain_bgn_id, interfaces.c.chain_end_id, unique=True)
Index('idx_interfaces_chain_end_id', interfaces.c.chain_end_id, interfaces.c.chain_bgn_id, unique=True)
Index('idx_interfaces_path', interfaces.c.path, postgresql_using='gist')

interface_comments = {
    "table": "Contains all binary protein-protein interactions.",
    "columns":
    {
        "interface_id": "Primary key of the interface.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "chain_bgn_id": "Primary key of the first chain.",
        "chain_end_id": "Primary key of the second chain.",
        "path": "ptree path in the form PDB/assembly serial number/I:PDB chain bgn ID-PDB chain end ID",
        "num_res_bgn": "Number of residues of the first chain that are part of this interface.",
        "num_res_end": "Number of residues of the second chain that are part of this interface.",
        "is_quaternary": "True if both chains are at identity, i.e. exist in the asymmetric unit.",
        "has_mod_res": "True if at least one modified residue if part of the interface.",
        "has_missing_atoms": "True if at least one residue has missing atoms."
    }
}

comment_on_table_elements(interfaces, interface_comments)

interface_residues = Table('interface_residues', metadata,
                           Column('interface_id', Integer, autoincrement=False),
                           Column('residue_bgn_id', Integer, autoincrement=False),
                           Column('residue_end_id', Integer, autoincrement=False),
                           schema=schema)

PrimaryKeyConstraint(interface_residues.c.interface_id, interface_residues.c.residue_bgn_id,
                     interface_residues.c.residue_end_id, deferrable=True, initially='deferred')
Index('idx_interface_residues_residue_bgn_id', interface_residues.c.residue_bgn_id)
Index('idx_interface_residues_residue_end_id', interface_residues.c.residue_end_id)


interface_residue_comments = {
    "table": "Contains the combinations of all residues in an interface that are interacting.",
    "columns":
    {
        "interface_id": "Primary key of the interface.",
        "residue_bgn_id": "Primary key of the interface residue of the first chain.",
        "residue_end_id": "Primary key of the interface residue of the second chain."
    }
}

comment_on_table_elements(interface_residues, interface_residue_comments)

#interface_atom_surface_areas = Table('interface_atom_surface_areas', metadata,
#                                     Column('interface_id', Integer, primary_key=True),
#                                     Column('atom_id', Integer, primary_key=True),
#                                     Column('asa_apo', Float(5,2)),
#                                     Column('asa_bound', Float(5,2)),
#                                     Column('asa_delta', Float(5,2)),
#                                     schema=schema)
#
#Index('idx_interface_atom_surface_areas_atom_id', interface_atom_surface_areas.c.atom_id)
#
#interface_atom_surface_area_comments = {
#    "table": "Solvent-accessible surface area changes for each atom in a protein-protein interface.",
#    "columns":
#    {
#        "interface_id": "Primary key of the interface.",
#        "atom_id": "Primary key of the atom that has a different solvent-exposed surface upon binding - can be a ligand or polymer atom.",
#        "asa_apo": "Solvent-accessible surface area of the atom in the apo state.",
#        "asa_bound": "Solvent-accessible surface area of the atom in the bound state.",
#        "asa_delta": "Change in solvent-accessible surface area between the apo and bound state."
#    }
#}
#
#comment_on_table_elements(interface_atom_surface_areas, interface_atom_surface_area_comments)
