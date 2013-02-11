"""
PDB residues are partitioned into residues, peptides, nucleotides and saccharides.
Residues belonging to ligands, water and residues without residue type remain in
the residues table.
"""

from sqlalchemy import (Boolean, CheckConstraint, Column, DDL, DefaultClause, Index,
                        Integer, String, Table, UniqueConstraint)
from sqlalchemy.event import listen
from sqlalchemy.schema import PrimaryKeyConstraint

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, comment_on_table_elements

# residues master table (won't be empty!)
residues = Table('residues', metadata,
                 Column('residue_id', Integer, nullable=False),
                 Column('biomolecule_id', Integer, nullable=False),
                 Column('chain_id', Integer, nullable=False),
                 Column('path', PTree, nullable=False),
                 Column('res_name', String(3), nullable=False),
                 Column('res_num', Integer, nullable=False),
                 Column('ins_code', String(1), nullable=False), # ' '
                 Column('entity_type_bm', Integer, nullable=False),
                 Column('is_disordered', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                 Column('is_incomplete', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                 schema=schema)

PrimaryKeyConstraint(residues.c.residue_id, deferrable=True, initially='deferred')
Index('idx_residues_biomolecule_id', residues.c.biomolecule_id)
Index('idx_residues_chain_id', residues.c.chain_id, residues.c.res_num, residues.c.ins_code)
Index('idx_residues_path', residues.c.path, postgresql_using='gist')

residue_comments = {
    "table": "Stores all residues in a PDB structure that are at least partially exposed to the surface and interact with another residue. This table is partitioned by entity type bit mask and this partition only has solvents and uknown residue types.",
    "columns":
    {
        "residue_id": "Primary key.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "chain_id": "Primary key of the parent chain.",
        "path": "ptree path in the form PDB/assembly serial number/PDB chain ID/PDB residue name[residue insertion code]`PDB residue number.",
        "res_name": "Three-letter name of the residue.",
        "res_num": "PDB residue number.",
        "ins_code": "PDB insertion code.",
        "entity_type_bm": "Entity type bitmask (the bits are solvent, ligand, saccharide, rna, dna, protein).",
        "is_disordered": "True if the residue has disordered atoms.",
        "is_incomplete": "True if the residue has missing atoms."
    }
}

comment_on_table_elements(residues, residue_comments)


### create the table partitions for the different polymer residue types


# polypeptide residues
peptides = Table('peptides', metadata,
                 Column('residue_id', Integer, nullable=False),
                 Column('biomolecule_id', Integer, nullable=False),
                 Column('chain_id', Integer, nullable=False),
                 Column('path', PTree, nullable=False),
                 Column('res_name', String(3), nullable=False),
                 Column('res_num', Integer, nullable=False),
                 Column('ins_code', String(1), nullable=False), # ' '
                 Column('entity_type_bm', Integer, nullable=False),
                 Column('is_disordered', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                 Column('is_incomplete', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                 Column('res_map_id', Integer),
                 Column('one_letter_code', String(1)),
                 Column('sstruct', String(1), DefaultClause('L'), nullable=False), # "'L'","'H'","'B'","'E'","'G'","'I'","'T'","'S'"
                 Column('cath', String(7)),
                 Column('px', Integer),
                 Column('is_non_std', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                 Column('is_modified', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                 Column('is_mutated', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                 CheckConstraint("entity_type_bm = 32 OR entity_type_bm = 34"),
                 schema=schema)

PrimaryKeyConstraint(peptides.c.residue_id, deferrable=True, initially='deferred')
Index('idx_peptides_biomolecule_id', peptides.c.biomolecule_id)
Index('idx_peptides_chain_id', peptides.c.chain_id, peptides.c.res_num, peptides.c.ins_code)
Index('idx_peptides_res_map_id', peptides.c.res_map_id) # CANNOT BE UNIQUE!
Index('idx_peptides_cath', peptides.c.cath)
Index('idx_peptides_px', peptides.c.px)
Index('idx_peptides_path', peptides.c.path, postgresql_using='gist')

# neccessary to drop tables with sqlalchemy
peptides.add_is_dependent_on(residues)

# add inheritance from master table through ddl
listen(peptides, "after_create",
       DDL("ALTER TABLE %(fullname)s INHERIT {schema}.residues".format(schema=schema)))

# DDL to create an insert rule on the master table
listen(metadata, "after_create",
       DDL("""
               CREATE OR REPLACE RULE peptide_insert_rule AS
                   ON INSERT TO {schema}.residues
                WHERE (entity_type_bm = 32 OR entity_type_bm = 34)
           DO INSTEAD INSERT INTO {schema}.peptides VALUES (NEW.*)
           """.format(schema=schema)))

# drop the rule on the master table because it depends on this table
listen(peptides, "before_drop",
       DDL("DROP RULE IF EXISTS {rule} ON {schema}.residues".format(schema=schema, rule='peptide_insert_rule')))

# add table comments
peptide_comments = {
    "table": "This partition stores all the residues that are part of a polypeptide sequence.",
    "columns":
    {
        "residue_id": "Primary key.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "chain_id": "Primary key of the parent chain.",
        "path": "ptree path in the form PDB/assembly serial number/PDB chain ID/PDB residue name[residue insertion code]`PDB residue number.",
        "res_name": "Three-letter name of the residue.",
        "res_num": "PDB residue number.",
        "ins_code": "PDB insertion code.",
        "entity_type_bm": "Entity type bitmask (the bits are solvent, ligand, saccharide, rna, dna, protein).",
        "is_disordered": "True if the residue has disordered atoms.",
        "is_incomplete": "True if the residue has missing atoms.",
        "res_map_id": "Primary key of this peptide in the SIFTS residue mapping.",
        "one_letter_code": "One-letter code of the amino acid. X is used in case of non-standard amino acids.",
        "sstruct": "Secondary structure DSSP code.",
        "cath": "CATH domain identifier.",
        "px": "SCOP px identifier.",
        "is_non_std": "True if the residue is not one of the 20 standard amino acids.",
        "is_modified": "True if this peptide is modified according to the SIFTS residue mapping.",
        "is_mutated": "True if this peptide is one of the 20 standard amino acids and the but differs from the canonical UniProt amino acid for this position."
    }
}

comment_on_table_elements(peptides, peptide_comments)

# dna/rna residues
nucleotides = Table('nucleotides', metadata,
                    Column('residue_id', Integer, nullable=False),
                    Column('biomolecule_id', Integer, nullable=False),
                    Column('chain_id', Integer, nullable=False),
                    Column('path', PTree, nullable=False),
                    Column('res_name', String(3), nullable=False),
                    Column('res_num', Integer, nullable=False),
                    Column('ins_code', String(1), nullable=False), # ' '
                    Column('entity_type_bm', Integer, nullable=False),
                    Column('is_disordered', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                    Column('is_incomplete', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                    CheckConstraint("entity_type_bm = 8 OR entity_type_bm = 16 OR entity_type_bm = 24"),
                    schema=schema)

PrimaryKeyConstraint(nucleotides.c.residue_id, deferrable=True, initially='deferred')
Index('idx_nucleotides_biomolecule_id', nucleotides.c.biomolecule_id)
Index('idx_nucleotides_chain_id', nucleotides.c.chain_id, nucleotides.c.res_num, nucleotides.c.ins_code)
Index('idx_nucleotides_path', nucleotides.c.path, postgresql_using='gist')

# neccessary to drop tables with sqlalchemy
nucleotides.add_is_dependent_on(residues)

# add inheritance from master table through ddl
listen(nucleotides, "after_create",
       DDL("ALTER TABLE %(fullname)s INHERIT {schema}.residues".format(schema=schema)))

# ddl to create an insert rule on the master table
listen(metadata, "after_create",
       DDL("""
               CREATE OR REPLACE RULE nucleotide_insert_rule AS
                   ON INSERT TO {schema}.residues
                WHERE (entity_type_bm = 8 OR entity_type_bm = 16 OR entity_type_bm = 24)
           DO INSTEAD INSERT INTO {schema}.nucleotides VALUES (NEW.*)
           """.format(schema=schema)))

# drop the rule on the master table because it depends on this table
listen(nucleotides, "before_drop",
       DDL("DROP RULE IF EXISTS {rule} ON {schema}.residues".format(schema=schema, rule='nucleotide_insert_rule')))

nucleotide_comments = {
    "table": "Stores all residues that are part of a oligonucleotide sequence.",
    "columns":
    {
        "residue_id": "Primary key.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "chain_id": "Primary key of the parent chain.",
        "path": "ptree path in the form PDB/assembly serial number/PDB chain ID/PDB residue name[residue insertion code]`PDB residue number.",
        "res_name": "Three-letter name of the residue.",
        "res_num": "PDB residue number.",
        "ins_code": "PDB insertion code.",
        "entity_type_bm": "Entity type bitmask (the bits are solvent, ligand, saccharide, rna, dna, protein).",
        "is_disordered": "True if the residue has disordered atoms.",
        "is_incomplete": "True if the residue has missing atoms."
    }
}

comment_on_table_elements(nucleotides, nucleotide_comments)

# saccharides
saccharides = Table('saccharides', metadata,
                    Column('residue_id', Integer, nullable=False),
                    Column('biomolecule_id', Integer, nullable=False),
                    Column('chain_id', Integer, nullable=False),
                    Column('path', PTree, nullable=False),
                    Column('res_name', String(3), nullable=False),
                    Column('res_num', Integer, nullable=False),
                    Column('ins_code', String(1), nullable=False), # ' '
                    Column('entity_type_bm', Integer, nullable=False),
                    Column('is_disordered', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                    Column('is_incomplete', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                    CheckConstraint("entity_type_bm = 4"),
                    schema=schema)

PrimaryKeyConstraint(saccharides.c.residue_id, deferrable=True, initially='deferred')
Index('idx_saccharides_biomolecule_id', saccharides.c.biomolecule_id)
Index('idx_saccharides_chain_id', saccharides.c.chain_id, saccharides.c.res_num, saccharides.c.ins_code)
Index('idx_saccharides_path', saccharides.c.path, postgresql_using='gist')

# NECCESSARY TO DROP TABLES WITH SQLALCHEMY
saccharides.add_is_dependent_on(residues)

# ADD INHERITANCE FROM MASTER TABLE THROUGH DDL
listen(saccharides, "after_create",
       DDL("ALTER TABLE %(fullname)s INHERIT {schema}.residues".format(schema=schema)))

# DDL TO CREATE AN INSERT RULE ON THE MASTER TABLE
listen(metadata, "after_create",
       DDL("""
               CREATE OR REPLACE RULE saccharide_insert_rule AS
                   ON INSERT TO {schema}.residues
                WHERE (entity_type_bm = 4)
           DO INSTEAD INSERT INTO {schema}.saccharides VALUES (NEW.*)
           """.format(schema=schema)))

# DROP THE RULE ON THE MASTER TABLE BECAUSE IT DEPENDS ON THIS TABLE
listen(saccharides, "before_drop",
       DDL("DROP RULE IF EXISTS {rule} ON {schema}.residues".format(schema=schema,
                                                          rule='saccharide_insert_rule')))

saccharide_comments = {
    "table": "Stores all residues that are part of a carbohydrate chain.",
    "columns":
    {
        "residue_id": "Primary key.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "chain_id": "Primary key of the parent chain.",
        "path": "ptree path in the form PDB/assembly serial number/PDB chain ID/PDB residue name[residue insertion code]`PDB residue number.",
        "res_name": "Three-letter name of the residue.",
        "res_num": "PDB residue number.",
        "ins_code": "PDB insertion code.",
        "entity_type_bm": "Entity type bitmask (the bits are solvent, ligand, saccharide, rna, dna, protein).",
        "is_disordered": "True if the residue has disordered atoms.",
        "is_incomplete": "True if the residue has missing atoms."
    }
}

comment_on_table_elements(saccharides, saccharide_comments)

#
residue_interaction_pairs = Table('residue_interaction_pairs', metadata,
                                  Column('biomolecule_id', Integer, nullable=False),
                                  Column('residue_bgn_id', Integer, nullable=False),
                                  Column('residue_end_id', Integer, nullable=False),
                                  Column('structural_interaction_type_bm', Integer, DefaultClause('0'), nullable=False),
                                  UniqueConstraint('residue_bgn_id', 'residue_end_id', name='residue_interaction_pairs_unique_interaction'),
                                  schema=schema)

Index('idx_residue_interaction_pairs_biomolecule_id', residue_interaction_pairs.c.biomolecule_id)
Index('idx_residue_interaction_pairs_residue_bgn_id', residue_interaction_pairs.c.residue_bgn_id)
Index('idx_residue_interaction_pairs_residue_end_id', residue_interaction_pairs.c.residue_end_id)

residue_interaction_pair_comments = {
    "table": "Contains all the residue pairs that are interacting with each other. This table is used as a shortcut to update interface/groove residues in particular.",
    "columns":
    {
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "residue_bgn_id": "Primary key of the residue.",
        "residue_end_id": "Primary key of the residue.",
        "structural_interaction_type": "Sum of the entity type bitmasks of the atom residues, e.g. 64 for protein-protein."
    }
}

comment_on_table_elements(residue_interaction_pairs, residue_interaction_pair_comments)
