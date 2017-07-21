from sqlalchemy import (Boolean, Column, DDL, DefaultClause, Float, Index, Integer,
                        String, Table, UniqueConstraint, CheckConstraint)
from sqlalchemy.schema import PrimaryKeyConstraint
from sqlalchemy.event import listen
from sqlalchemy.dialects.postgresql import ARRAY, array
from sqlalchemy.sql.elements import quoted_name 

from credovi import app
from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, Vector3D, comment_on_table_elements, create_partition_insert_trigger

# table to hold non-ring pi_groups
pi_groups = Table('pi_groups',metadata,
                   Column('pi_id', Integer, nullable=False),
                   Column('biomolecule_id', Integer, nullable=False),
                   #Column('residue_ids', ARRAY(Integer), nullable=False),
                   #Column('path', ARRAY(PTree)),
                   Column('pi_serial', Integer, nullable=False),
                   #Column('pi_number', Integer),
                   Column('size', Integer, nullable=False),
                   Column('centroid', Vector3D),
                   Column('normal', Vector3D),
                   #Column('is_hetero_aromatic', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                   schema=schema)

PrimaryKeyConstraint(pi_groups.c.pi_id, deferrable=True, initially='deferred')
Index('idx_pi_biomolecule_id', pi_groups.c.biomolecule_id, pi_groups.c.pi_serial, unique=True)
#Index('idx_pi_path', pi_groups.c.path, postgresql_using='gist')
#Index('idx_pi_resid', pi_groups.c.residue_ids, postgresql_using='gin')

pi_groups_comments = {
    "table": "Contains the non-ring pi groups found in structures.",
    "columns":
    {
        "pi_id": "Primary Key of the pi system.",
        "biomolecule_id": "Foreign key of the parent biomolecule.",
        "pi_serial": "Pi system serial number inside the parent biomolecule.",
        #"pi_number": "Number of the pi system inside the residue it is part of.",
        "size": "Number of atoms in the pi group.",
        "centroid": "The Centroid of the pi group as vector3d type.",
        "normal": "The normal of the pi group as vector3d type.",
        #"is_hetero_aromatic": "True if the aromatic ring contains atoms other than carbon."
    }
}

comment_on_table_elements(pi_groups, pi_groups_comments)


## table to hold pi group atoms
CURRENT_BIOMOL_MAX      = app.config.get('schema','current_biomol_max')
PI_GROUP_PARTITION_SIZE = app.config.get('schema','pi_group_partition_size')  # In terms of biomolecules

PI_GROUP_INS_RULE_DDL = """
                        CREATE OR REPLACE RULE {rule} AS
                            ON INSERT TO {schema}.{master_table}
                        WHERE (biomolecule_id > {part_bound_low} AND biomolecule_id <= {part_bound_high})
                        DO INSTEAD INSERT INTO {schema}.{table} VALUES (NEW.*)
                        """

## pi group to residues map
pi_group_res = Table('pi_group_residues', metadata,
                     Column('pi_id', Integer, nullable=False, autoincrement=False),
                     Column('biomolecule_id', Integer, nullable=False),
                     Column('residue_id', Integer, nullable=False),
                     Column('path', PTree, nullable=False),
                     schema=schema)

PrimaryKeyConstraint(pi_group_res.c.pi_id, pi_group_res.c.residue_id, deferrable=True, initially='deferred')
#
#Index('idx_pi_res_path', pi_group_res.c.path, postgresql_using='gist')

pi_group_res_comments = {
    "table": "Maps pi groups residues to their ids and paths",
    "columns":
    {
        "pi_id": "Primary Key of the pi system.",
        "biomolecule_id": "Foreign key of the parent biomolecule.",
        "residue_id": "Foreign key of the residue(s) the pi group is part of.",
        "path": "ptree path in the form PDB/assembly serial number/PDB chain ID/PDB residue name[residue insertion code]`PDB residue number/PI:pi number.",
    }
}

comment_on_table_elements(pi_group_res, pi_group_res_comments)
create_partition_insert_trigger(pi_group_res, PI_GROUP_PARTITION_SIZE)


pi_atoms = Table('pi_group_atoms', metadata,
                  Column('pi_atom_id', Integer, primary_key=True),
                  Column('pi_id', Integer, nullable=False),
                  Column('biomolecule_id', Integer, nullable=False),
                  Column('atom_id', Integer, nullable=False),
                  schema=schema)

pi_atom_comments = {
    "table": "Contains the atoms that a given pi system from the pi_groups table consists of.",
    "columns":
    {
        "pi_atom_id": "Primary key of the pi systems atom.",
        "pi_id": "Primary key of the pi system this atom is part of.",
        "biomolecule_id": "Foreign key of the parent biomolecule.",
        "atom_id": "Primary key of the atom."
    }
}

comment_on_table_elements(pi_atoms, pi_atom_comments)
create_partition_insert_trigger(pi_atoms, PI_GROUP_PARTITION_SIZE)


# interactions between pi groups (and rings)
pi_interactions = Table('pi_interactions', metadata,
                        Column('pi_interaction_id', Integer, primary_key=True),
                        Column('biomolecule_id', Integer, nullable=False),
                        Column('pi_bgn_id', Integer, nullable=False),
                        Column('pi_bgn_is_ring', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                        Column('pi_end_id', Integer, nullable=False),
                        Column('pi_end_is_ring', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                        # Column('closest_atom_bgn_id', Integer),
                        # Column('closest_atom_end_id', Integer),
                        Column('distance', Float(3,2), nullable=False),
                        # Column('closest_atom_distance', Float(3,2)),
                        Column('dihedral', Float(5,2), nullable=False),
                        Column('theta', Float(5,2), nullable=False),
                        Column('iota', Float(5,2), nullable=False),
                        Column('interaction_type', String(2)), # "'FF'", "'OF'", "'EE'", "'FT'", "'OT'", "'ET'", "'FE'", "'OE'", "'EF'")
                        schema=schema)

pi_interaction_comments = {
    "table": "Contains the parameters of the interaction between two pi groups.",
    "columns":
    {
        "pi_interaction_id": "Primary key of the pi interaction.",
        "biomolecule_id": "Foreign key of the parent biomolecule.",
        "pi_bgn_id": "Primary key of the first pi group.",
        "pi_bgn_is_ring": "Whether the first pi group is an aromatic ring or a different type.",
        "pi_end_id": "Primary key of the second pi group.",
        "pi_end_is_ring": "Whether the second pi group is an aromatic ring or a different type.",
        # "closest_atom_bgn_id": "Primary key of the closest atom of the first pi group.",
        # "closest_atom_end_id": "Primary key of the closest atom of the second pi group.",
        "distance": "Distance between the pi group centroids.",
        # "closest_atom_distance": "Distance between the closest atoms of the pi groups.",
        "dihedral": "Dihedral angle between the two planes (normals).",
        "theta": "Signed angle in degrees between the normal of the first pi group and the vector between the two centroids.",
        "iota": "Signed angle in degrees between the normal of the second pi group and the vector between the two centroids.",
        "interaction_type": "Classification of the pi interaction geometry."
    }
}

# Index('idx_pi_interactions_biomol_id_atom_bgn_id', pi_interactions.c.biomolecule_id, pi_interactions.c.closest_atom_bgn_id)
# Index('idx_pi_interactions_biomol_id_atom_end_id', pi_interactions.c.biomolecule_id, pi_interactions.c.closest_atom_end_id)
# Index('idx_pi_interactions_pi_bgn_id', pi_interactions.c.pi_bgn_is_ring, pi_interactions.c.pi_bgn_id)
# Index('idx_pi_interactions_pi_end_id', pi_interactions.c.pi_end_is_ring, pi_interactions.c.pi_end_id)
comment_on_table_elements(pi_interactions, pi_interaction_comments)
create_partition_insert_trigger(pi_interactions, PI_GROUP_PARTITION_SIZE)



# create new partitions for every X biomolecules
partitions = range(0, CURRENT_BIOMOL_MAX+PI_GROUP_PARTITION_SIZE, PI_GROUP_PARTITION_SIZE)


for part_bound_low, part_bound_high in zip(partitions[:-1], partitions[1:]):
    ## RESIDUES
    residues_tablename = 'pi_group_residues_biomol_le_{0}'.format(part_bound_high)
    # residues_rulename = residues_tablename + '_insert'  ## TRIGGERS, NOT RULES!

    residues_partition = Table(residues_tablename, metadata,
                               Column('pi_id', Integer, nullable=False, autoincrement=False),
                               Column('biomolecule_id', Integer, nullable=False),
                               Column('residue_id', Integer, nullable=False),
                               Column('path', PTree, nullable=False),
                               CheckConstraint("biomolecule_id > {0} AND biomolecule_id <= {1}".format(part_bound_low, part_bound_high)),
                               postgresql_inherits=quoted_name(pi_group_res.fullname, False),  # new SQLAlchemy 1.0 feature
                               schema=schema)
                               
    PrimaryKeyConstraint(residues_partition.c.pi_id, residues_partition.c.residue_id, deferrable=True, initially='deferred')
    Index('idx_{}_path'.format(residues_tablename), residues_partition.c.path, postgresql_using='gist')

    # neccessary to drop tables with sqlalchemy
    residues_partition.add_is_dependent_on(pi_group_res)

    # # add inheritance from master table through ddl - NO LONGER NECESSARY UNDER SQLALCHEMY 1.0
    listen(residues_partition, "after_create",
           DDL("ALTER TABLE %(fullname)s INHERIT {schema}.pi_group_residues".format(schema=schema)))

    ## ddl to create an insert rule on the master table  # DEPRECATED
    # listen(metadata, "after_create",
           # DDL(PI_GROUP_INS_RULE_DDL.format(schema=schema, master_table='pi_group_residues',
           #                                  table=residues_tablename,
           #                                  rule=residues_rulename, part_bound_low=part_bound_low,
           #                                  part_bound_high=part_bound_high))
           # )

    # # drop the rules on the master table because they depend on the partitions # DEPRECATED
    # listen(residues_partition, "before_drop",
    #        DDL("DROP RULE IF EXISTS {rule} ON {schema}.pi_group_residues".format(schema=schema, rule=residues_rulename)))

    comment_on_table_elements(residues_partition, pi_group_res_comments)


    ## ATOMS
    atoms_tablename = 'pi_group_atoms_biomol_le_{0}'.format(part_bound_high)
    #atoms_rulename = atoms_tablename + '_insert'  #  TRIGGERS, NOT RULES!

    atoms_partition = Table(atoms_tablename, metadata,
                      Column('pi_atom_id', Integer, primary_key=True),
                      Column('pi_id', Integer, nullable=False, index=True),
                      Column('biomolecule_id', Integer, nullable=False, index=True),
                      Column('atom_id', Integer, nullable=False, index=True),
                      CheckConstraint("biomolecule_id > {0} AND biomolecule_id <= {1}".format(part_bound_low, part_bound_high)),
                      postgresql_inherits=quoted_name(pi_atoms.fullname, False),  # new SQLAlchemy 1.0 feature 
                      schema=schema)

    # neccessary to drop tables with sqlalchemy
    atoms_partition.add_is_dependent_on(pi_atoms)

    # # add inheritance from master table through ddl # DEPRECATED UNDER SQLALCHEMY 1.0+
    # listen(atoms_partition, "after_create",
            # DDL("ALTER TABLE %(fullname)s INHERIT {schema}.pi_group_atoms".format(schema=schema)))

    ## ddl to create an insert rule on the master table # DEPRECATED
    # listen(metadata, "after_create",
           # DDL(PI_GROUP_INS_RULE_DDL.format(schema=schema, master_table='pi_group_atoms', table=atoms_tablename,
                                            # rule=atoms_rulename, part_bound_low=part_bound_low,
                                            # part_bound_high=part_bound_high)))
    
    # # drop the rules on the master table because they depend on the partitions # DEPRECATED
    # listen(atoms_partition, "before_drop",
           # DDL("DROP RULE IF EXISTS {rule} ON {schema}.pi_group_atoms".format(schema=schema, rule=atoms_rulename)))

    comment_on_table_elements(atoms_partition, pi_atom_comments)


    ### PI Interactions
    pi_inter_tablename = 'pi_interactions_biomol_le_{0}'.format(part_bound_high)

    pi_inter_partition = Table(pi_inter_tablename, metadata,
                                Column('pi_interaction_id', Integer, primary_key=True),
                                Column('biomolecule_id', Integer, index=True, nullable=False),
                                Column('pi_bgn_id', Integer, nullable=False),
                                Column('pi_bgn_is_ring', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                                Column('pi_end_id', Integer, nullable=False),
                                Column('pi_end_is_ring', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                                # Column('closest_atom_bgn_id', Integer),
                                # Column('closest_atom_end_id', Integer),
                                Column('distance', Float(3,2), nullable=False),
                                # Column('closest_atom_distance', Float(3,2)),
                                Column('dihedral', Float(5,2), nullable=False),
                                Column('theta', Float(5,2), nullable=False),
                                Column('iota', Float(5,2), nullable=False),
                                Column('interaction_type', String(2)), # "'FF'", "'OF'", "'EE'", "'FT'", "'OT'", "'ET'", "'FE'", "'OE'", "'EF'")
                                CheckConstraint("biomolecule_id > {0} AND biomolecule_id <= {1}".format(part_bound_low, part_bound_high)),
                                postgresql_inherits=quoted_name(pi_atoms.fullname, False),  # new SQLAlchemy 1.0 feature
                                schema=schema)

    pi_inter_partition.add_is_dependent_on(pi_interactions)
    
    ## DEPRECATED UNDER SQLALCHEMY 1.0+
    # listen(pi_inter_partition, "after_create",
            # DDL("ALTER TABLE %(fullname)s INHERIT {schema}.pi_interactions".format(schema=schema)))

    comment_on_table_elements(pi_inter_partition, pi_interaction_comments)
