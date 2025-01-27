"""
The atoms table is partitioned by biomolecule_id. Every partition will hold the
atoms of 10,000 biomolecules.
"""
from sqlalchemy import (Boolean, CheckConstraint, Column, DDL, Float, Index, Integer,
                        SmallInteger, String, Table, DefaultClause)
from sqlalchemy.event import listen
from sqlalchemy.sql.elements import quoted_name
from sqlalchemy.schema import PrimaryKeyConstraint

from credovi import app
from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, Vector3D, comment_on_table_elements, create_partition_insert_trigger

CURRENT_BIOMOL_MAX   = app.config.get('schema','current_biomol_max')
ATOMS_PARTITION_SIZE = app.config.get('schema','atoms_partition_size')  # In terms of biomolecules

ATOMS_INS_RULE_DDL = """
                        CREATE OR REPLACE RULE {rule} AS
                            ON INSERT TO {schema}.atoms
                         WHERE (biomolecule_id > {part_bound_low} AND biomolecule_id <= {part_bound_high})
                    DO INSTEAD INSERT INTO {schema}.{table} VALUES (NEW.*)
                    """

# master table
atoms = Table('atoms', metadata,
              Column('atom_id', Integer, primary_key=True),
              Column('biomolecule_id', Integer, nullable=False),
              Column('residue_id', Integer, nullable=False),
              Column('path', PTree),
              Column('atom_serial', Integer, nullable=False),
              Column('group_pdb', String(7), nullable=False), # HETATM/ATOM
              Column('atom_name', String(4), nullable=False),
              Column('alt_loc', String(1), nullable=False),
              Column('coords', Vector3D, nullable=False),
              Column('occupancy', Float(3,2), nullable=False),
              Column('b_factor', Float(4,2), nullable=False),
              Column('element', String(2)),
              Column('hyb', SmallInteger, nullable=False),
              Column('tripos_atom_type', String(5)),
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
              schema=schema)

comments = {
    "table": "Stores ALL PDB atoms belonging that are part of a residue.",
    "columns":
    {
        "atom_id": "Primary key of the atom.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "residue_id": "Primary key of the parent residue.",
        "path": "PyMOL-like selection syntax for this atom.",
        "atom_serial": "PDB serial number of the atom.",
        "group_pdb": "Either ATOM or HETATM.",
        "atom_name": "PDB Atom name.",
        "alt_loc": "Alternate location code if disordered otherwise single space.",
        "coords": "Coordinates of the atom as vector3d type.",
        "occupancy": "Atom occupancy.",
        "b_factor": "Atom B-factor.",
        "element": "Element of the atom",
        "hyb": "OEChem atom hybridisation.",
        "tripos_atom_type": "Tripos atom type according to OEChem.",
        "is_donor": "",
        "is_acceptor": "",
        "is_aromatic": "",
        "is_weak_acceptor": "",
        "is_weak_donor": "",
        "is_hydrophobe": "",
        "is_metal": "",
        "is_pos_ionisable": "",
        "is_neg_ionisable": "",
        "is_xbond_donor": "",
        "is_xbond_acceptor": "",
        "is_carbonyl_oxygen": "",
        "is_carbonyl_carbon": ""
    }
}

comment_on_table_elements(atoms, comments)
create_partition_insert_trigger(atoms, ATOMS_PARTITION_SIZE)

# create new paritions for every X biomolecules
partitions = range(0, CURRENT_BIOMOL_MAX+ATOMS_PARTITION_SIZE, ATOMS_PARTITION_SIZE)

for part_bound_low, part_bound_high in zip(partitions[:-1], partitions[1:]):
    tablename = 'atoms_biomol_le_{0}'.format(part_bound_high)
    rulename = tablename + '_insert'

    partition = Table(tablename, metadata,
                      Column('atom_id', Integer, primary_key=True),
                      Column('biomolecule_id', Integer, nullable=False),
                      Column('residue_id', Integer, nullable=False),
                      Column('path', PTree),
                      Column('atom_serial', Integer, nullable=False),
                      Column('group_pdb', String(7), nullable=False), # HETATM/ATOM
                      Column('atom_name', String(4), nullable=False),
                      Column('alt_loc', String(1), nullable=False),
                      Column('coords', Vector3D, nullable=False),
                      Column('occupancy', Float(3,2), nullable=False),
                      Column('b_factor', Float(4,2), nullable=False),
                      Column('element', String(2)),
                      Column('hyb', SmallInteger, nullable=False),
                      Column('tripos_atom_type', String(5)),
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
                      CheckConstraint("biomolecule_id > {0} AND biomolecule_id <= {1}".format(part_bound_low, part_bound_high)),
                      postgresql_inherits=quoted_name(atoms.fullname, False),  # new SQLAlchemy 1.0 feature
                      schema=schema)

    Index('idx_{0}_atom'.format(tablename), partition.c.residue_id, partition.c.atom_name,
          partition.c.alt_loc, unique=True)

    Index('idx_{0}_serial'.format(tablename), partition.c.biomolecule_id,
          partition.c.atom_serial, unique=True)

    # neccessary to drop tables with sqlalchemy
    partition.add_is_dependent_on(atoms)

    # # add inheritance from master table through ddl
    # listen(partition, "after_create",
    #        DDL("ALTER TABLE %(fullname)s INHERIT {schema}.atoms".format(schema=schema)))
    #
    # # ddl to create an insert rule on the master table
    # listen(metadata, "after_create",
    #        DDL(ATOMS_INS_RULE_DDL.format(schema=schema, table=tablename,
    #                                      rule=rulename, part_bound_low=part_bound_low,
    #                                      part_bound_high=part_bound_high)))
    #
    # # drop the rules on the master table because they depend on the partitions
    # listen(partition, "before_drop",
    #        DDL("DROP RULE IF EXISTS {rule} ON {schema}.atoms".format(schema=schema, rule=rulename)))

    comment_on_table_elements(partition, comments)