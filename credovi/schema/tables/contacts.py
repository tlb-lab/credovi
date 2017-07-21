"""
The contacts table is too large and has to be partitioned into more sizeable chunks.
Currently, the table is partitioned using the biomolecule_id column - every partition
will hold 5000 biomolecules.
"""

from sqlalchemy import Boolean, CheckConstraint, Column, DDL, Float, Index, Integer, Table, DefaultClause
from sqlalchemy.event import listen
from sqlalchemy.schema import PrimaryKeyConstraint

from credovi import app
from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import comment_on_table_elements

CURRENT_BIOMOL_MAX      = app.config.get('schema','current_biomol_max')
CONTACTS_PARTITION_SIZE = app.config.get('schema','contacts_partition_size')

CONTACTS_INS_RULE_DDL   = """
                            CREATE OR REPLACE RULE {rule} AS
                                ON INSERT TO {schema}.contacts
                             WHERE (biomolecule_id > {part_bound_low} AND biomolecule_id <= {part_bound_high})
                        DO INSTEAD INSERT INTO {schema}.{table} VALUES (NEW.*)
                        """

# master table
contacts = Table('contacts', metadata,
                 Column('contact_id', Integer, primary_key=True, nullable=False),
                 Column('biomolecule_id', Integer, nullable=False),
                 Column('atom_bgn_id', Integer, nullable=False),
                 Column('atom_end_id', Integer, nullable=False),
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
                 schema=schema)

comments = {
    "table": "Contains the contacts between all interacting atoms.",
    "columns":
    {
        "contact_id": "Primary key of the contact.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "atom_bgn_id": "Primary key of the first atom.",
        "atom_end_id": "Primary key of the second atom.",
        "distance": "Distance between the atoms.",
        "structural_interaction_type_bm": "Bit mask created by concatenating the entity type bit masks of the interacting atoms.",
        "is_clash": "",
        "is_same_entity":"True if both atoms belong to the same entity.",
        "is_covalent": "",
        "is_vdw_clash": "",
        "is_vdw": "",
        "is_proximal": "",
        "is_hbond": "",
        "is_weak_hbond": "",
        "is_xbond": "",
        "is_ionic": "",
        "is_metal_complex": "",
        "is_aromatic": "",
        "is_hydrophobic": "",
        "is_carbonyl": ""
        }
}

comment_on_table_elements(contacts, comments)

# create new paritions for every 5000 biomolecules
partitions = range(0, CURRENT_BIOMOL_MAX+CONTACTS_PARTITION_SIZE, CONTACTS_PARTITION_SIZE)

for part_bound_low, part_bound_high in zip(partitions[:-1], partitions[1:]):
    tablename = 'contacts_biomol_le_{0}'.format(part_bound_high)
    rulename = tablename + '_insert'

    partition = Table(tablename, metadata,
                    Column('contact_id', Integer, primary_key=True, nullable=False),
                    Column('biomolecule_id', Integer, nullable=False),
                    Column('atom_bgn_id', Integer, nullable=False),
                    Column('atom_end_id', Integer, nullable=False),
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
                      CheckConstraint("biomolecule_id > {0} AND biomolecule_id <= {1}".format(part_bound_low, part_bound_high)),
                      schema=schema)

    Index('idx_{0}_biomolecule_id'.format(tablename), partition.c.biomolecule_id)
    #Index('idx_{0}_atom_bgn_id'.format(tablename), partition.c.atom_bgn_id)
    #Index('idx_{0}_atom_end_id'.format(tablename), partition.c.atom_end_id)

    # neccessary to drop tables with sqlalchemy
    partition.add_is_dependent_on(contacts)

    # add inheritance from master table through DDL
    listen(partition, "after_create",
           DDL("ALTER TABLE %(fullname)s INHERIT {schema}.contacts".format(schema=schema)))

    # DDL to create an insert rule on the master table
    listen(metadata, "after_create",
           DDL(CONTACTS_INS_RULE_DDL.format(schema=schema, table=tablename,
                                            rule=rulename, part_bound_low=part_bound_low,
                                            part_bound_high=part_bound_high)))

    # drop the rules on the master table because they depend on the partitions
    listen(partition, "before_drop",
       DDL("DROP RULE IF EXISTS {rule} ON {schema}.contacts".format(schema=schema, rule=rulename)))

    comment_on_table_elements(partition, comments)
