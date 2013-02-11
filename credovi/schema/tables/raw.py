"""
This script contains all the tables that are used for loading the data dumps created
for each PDB file.
"""

from sqlalchemy import Boolean, Column, DefaultClause, Float, Index, Integer, String, Table, Text
from sqlalchemy.dialects.postgresql import ARRAY, DOUBLE_PRECISION

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import Vector3D

# RAW CHAINS
raw_chains = Table('raw_chains', metadata,
                   Column('pdb', String(4), nullable=False, primary_key=True),
                   Column('assembly_serial', Integer, nullable=False, primary_key=True, autoincrement=False),
                   Column('entity_serial', Integer, nullable=False, autoincrement=False, primary_key=True),
                   Column('pdb_chain_id', String(1), nullable=False),
                   Column('pdb_chain_asu_id', String(1), nullable=False),
                   Column('chain_type', String(50)),
                   Column('rotation', ARRAY(DOUBLE_PRECISION)),
                   Column('translation', ARRAY(DOUBLE_PRECISION)),
                   Column('is_at_identity', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                   schema=schema, prefixes=['unlogged'])

# RAW LIGANDS
raw_ligands = Table('raw_ligands', metadata,
                    Column('pdb', String(4), nullable=False, primary_key=True),
                    Column('assembly_serial', Integer, nullable=False, primary_key=True, autoincrement=False),
                    Column('entity_serial', Integer, nullable=False, autoincrement=False, primary_key=True),
                    Column('pdb_chain_id', String(1), nullable=False),
                    Column('res_num', Integer),
                    Column('ligand_name', String(64), nullable=False),
                    Column('num_hvy_atoms', Integer),
                    #Column('ism', Text),
                    schema=schema, prefixes=['unlogged'])

# RAW ATOM TABLE
raw_atoms = Table('raw_atoms', metadata,
                  Column('pdb', String(4), nullable=False, primary_key=True),
                  Column('assembly_serial', Integer, nullable=False, primary_key=True, autoincrement=False),
                  Column('atom_serial', Integer, nullable=False, primary_key=True, autoincrement=False),
                  Column('group_pdb', String(7), nullable=False),
                  Column('atom_name', String(4), nullable=False),
                  Column('alt_loc', String(1), nullable=False),
                  Column('res_name', String(3), nullable=False),
                  Column('pdb_chain_id', String(1), nullable=False),
                  Column('pdb_chain_asu_id', String(1), nullable=False),
                  Column('res_num', Integer, nullable=False),
                  Column('ins_code', String(1), nullable=False),
                  Column('coords', Vector3D, nullable=False),
                  Column('occupancy', Float(6,2), nullable=False),
                  Column('b_factor', Float(6,2), nullable=False),
                  Column('element', String(2)),
                  Column('hyb', Integer, nullable=False),
                  Column('tripos_atom_type', String(5)),
                  Column('entity_serial', Integer, nullable=False),
                  Column('entity_type_bm', Integer, nullable=False),
                  Column('is_donor', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_acceptor', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_aromatic', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_weak_acceptor', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_weak_donor', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_hydrophobe', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_metal', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_pos_ionisable', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_neg_ionisable', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_xbond_donor', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_xbond_acceptor', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_carbonyl_oxygen', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  Column('is_carbonyl_carbon', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                  schema=schema, prefixes=['unlogged'])

Index('idx_raw_atoms', raw_atoms.c.pdb, raw_atoms.c.assembly_serial,
      raw_atoms.c.pdb_chain_id, raw_atoms.c.res_num, raw_atoms.c.ins_code,
      raw_atoms.c.atom_name)

raw_contacts = Table('raw_contacts', metadata,
                     Column('pdb', String(4), primary_key=True),
                     Column('assembly_serial', Integer, primary_key=True, autoincrement=False),
                     Column('atom_bgn_serial', Integer, primary_key=True, autoincrement=False),
                     Column('atom_end_serial', Integer, primary_key=True, autoincrement=False),
                     Column('distance', Float(3,2), nullable=False),
                     Column('structural_interaction_type_bm', Integer, DefaultClause('0'), nullable=False),
                     Column('is_same_entity', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_clash', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_covalent', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_vdw_clash', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_vdw', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_proximal', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_hbond', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_weak_hbond', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_xbond', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_ionic', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_metal_complex', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_aromatic', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_hydrophobic', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_carbonyl', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     schema=schema, prefixes=['unlogged'])

raw_rings = Table('raw_aromatic_rings', metadata,
                  Column('pdb', String(4), nullable=False, primary_key=True),
                  Column('assembly_serial', Integer, nullable=False, primary_key=True, autoincrement=False),
                  Column('ring_serial', Integer, nullable=False, primary_key=True, autoincrement=False),
                  Column('atom_serial', Integer, nullable=False, primary_key=True, autoincrement=False),
                  schema=schema, prefixes=['unlogged'])

raw_binding_site_atom_surface_areas = Table('raw_binding_site_atom_surface_areas', metadata,
                 Column('pdb', String(4), nullable=False, primary_key=True),
                 Column('assembly_serial', Integer, nullable=False, autoincrement=False, primary_key=True),
                 Column('entity_serial', Integer, nullable=False, autoincrement=False, primary_key=True), # USED TO JOIN THE LIGANDS TABLE
                 Column('atom_serial', Integer, nullable=False, autoincrement=False, primary_key=True), # USED TO JOIN THE ATOMS
                 Column('asa_apo', Float(6,3)),
                 Column('asa_bound', Float(6,3)),
                 Column('asa_delta', Float(6,3)),
                 schema=schema, prefixes=['unlogged'])

Index('idx_raw_binding_site_atom_surface_areas_entity_serial',
      raw_binding_site_atom_surface_areas.c.pdb,
      raw_binding_site_atom_surface_areas.c.assembly_serial,
      raw_binding_site_atom_surface_areas.c.entity_serial)

Index('idx_raw_binding_site_atom_surface_areas_atom_serial',
      raw_binding_site_atom_surface_areas.c.pdb,
      raw_binding_site_atom_surface_areas.c.assembly_serial,
      raw_binding_site_atom_surface_areas.c.atom_serial)

raw_interface_atom_surface_areas = Table('raw_interface_atom_surface_areas', metadata,
             Column('pdb', String(4), nullable=False),
             Column('assembly_serial', Integer, nullable=False, autoincrement=False),
             Column('pdb_chain_bgn_id', String(1), nullable=False),
             Column('pdb_chain_end_id', String(1), nullable=False),
             Column('pdb_chain_id', String(1), nullable=False),
             Column('atom_serial', Integer, nullable=False, autoincrement=False),
             Column('asa_apo', Float(5,2)),
             Column('asa_bound', Float(5,2)),
             Column('asa_delta', Float(5,2)),
             schema=schema)

Index('idx_raw_interface_atom_surface_areas_interface',
      raw_interface_atom_surface_areas.c.pdb,
      raw_interface_atom_surface_areas.c.assembly_serial,
      raw_interface_atom_surface_areas.c.pdb_chain_bgn_id,
      raw_interface_atom_surface_areas.c.pdb_chain_end_id)

Index('idx_raw_interface_atom_surface_areas_pdb_chain_id',
      raw_interface_atom_surface_areas.c.pdb,
      raw_interface_atom_surface_areas.c.assembly_serial,
      raw_interface_atom_surface_areas.c.pdb_chain_id,
      raw_interface_atom_surface_areas.c.atom_serial)
