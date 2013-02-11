from sqlalchemy import Column, DefaultClause, Index, Integer, String, Table, DefaultClause
from sqlalchemy.schema import PrimaryKeyConstraint

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, comment_on_table_elements

biomolecules = Table('biomolecules', metadata,
                     Column('biomolecule_id', Integer, nullable=False),
                     Column('structure_id', Integer, nullable=False),
                     Column('path', PTree, nullable=False),
                     Column('assembly_serial', Integer, nullable=False),
                     Column('assembly_type', String(32)), # e.g. monomeric
                     Column('conformational_state_bm', Integer, DefaultClause('0'), nullable=False),
                     Column('structural_interaction_bm', Integer, DefaultClause('0'), nullable=False),
                     Column('num_chains', Integer, DefaultClause('0'), nullable=False),
                     Column('num_ligands', Integer, DefaultClause('0'), nullable=False),
                     Column('num_atoms', Integer, DefaultClause('0'), nullable=False),
                     schema=schema)

PrimaryKeyConstraint(biomolecules.c.biomolecule_id, deferrable=True, initially='deferred')
Index('idx_biomolecules_structures_id', biomolecules.c.structure_id, biomolecules.c.assembly_serial, unique=True)
Index('idx_biomolecules_path', biomolecules.c.path, postgresql_using='gist')

comments = {
    "table": "Stores the biological assemblies that are generated from an asymmetric unit.",
    "columns":
    {
        "biomolecule_id": "Primary key.",
        "structure_id": "Primary key of the parent structure.",
        "path": "ptree path in the form PDB/assembly serial number.",
        "assembly_serial": "Serial number of the assembly. Can be > 1 if asymmetric unit contains more than one biological assembly.",
        "assembly_type": "Assembly type of the biomolecule, e.g. Monomeric.",
        "conformational_state_bm": "Bit mask of the conformational state of the biomolecule.",
        "structural_interaction_bm": "Bit mask of the kinds of structural interactions that are present.",
        "num_chains": "Number of chains in biological assembly.",
        "num_ligands": "Number of ligands in biological assembly.",
        "num_atoms": "Number of atoms in biological assembly."
    }
}

comment_on_table_elements(biomolecules, comments)
