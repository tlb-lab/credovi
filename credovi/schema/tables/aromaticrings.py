from sqlalchemy import Boolean, Column, DDL, DefaultClause, Float, Index, Integer, String, Table, UniqueConstraint

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, Vector3D, comment_on_table_elements

# table to hold aromatic rings
aromatic_rings = Table('aromatic_rings',metadata,
                       Column('aromatic_ring_id', Integer, primary_key=True),
                       Column('biomolecule_id', Integer, nullable=False),
                       Column('residue_id', Integer, nullable=False),
                       Column('path', PTree),
                       Column('ring_serial', Integer, nullable=False),
                       Column('ring_number', Integer),
                       Column('size', Integer, nullable=False),
                       Column('centroid', Vector3D),
                       Column('normal', Vector3D),
                       Column('is_hetero_aromatic', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                       schema=schema)

Index('idx_aromatic_rings', aromatic_rings.c.residue_id, aromatic_rings.c.ring_serial, unique=True)
Index('idx_aromatic_rings_biomolecule_id', aromatic_rings.c.biomolecule_id, aromatic_rings.c.ring_serial, unique=True)
Index('idx_aromatic_rings_path', aromatic_rings.c.path, postgresql_using='gist')

aromatic_ring_comments = {
    "table": "Contains the aromatic rings found in aromatic residues from the residues table.",
    "columns":
    {
        "aromatic_ring_id": "Primary Key of the aromatic ring.",
        "biomolecule_id": "Foreign key of the parent biomolecule.",
        "residue_id": "Foreign key of the residue the aromatic ring is part of.",
        "path": "ptree path in the form PDB/assembly serial number/PDB chain ID/PDB residue name[residue insertion code]`PDB residue number/AR:ring number.",
        "ring_serial": "Aromatic ring serial number inside the parent biomolecule.",
        "ring_number": "Number of the aromatic ring inside the residue it is part of.",
        "size": "Number of atoms in the aromatic ring.",
        "centroid": "The Centroid of the aromatic ring as vector3d type.",
        "normal": "The normal of the aromatic ring as vector3d type.",
        "is_hetero_aromatic": "True if the aromatic ring contains atoms other than carbon."
    }
}

comment_on_table_elements(aromatic_rings, aromatic_ring_comments)

aromatic_ring_atoms = Table('aromatic_ring_atoms', metadata,
                            Column('aromatic_ring_atom_id', Integer, primary_key=True),
                            Column('aromatic_ring_id', Integer, nullable=False),
                            Column('atom_id', Integer, nullable=False),
                            schema=schema)

Index('idx_aromatic_ring_atoms_aromatic_ring_id', aromatic_ring_atoms.c.aromatic_ring_id)
Index('idx_aromatic_ring_atoms_atom_id', aromatic_ring_atoms.c.atom_id)

aromatic_ring_atom_comments = {
    "table": "Contains the atoms that a given aromatic ring from the aromatic_rings table consists of.",
    "columns":
    {
        "aromatic_ring_atom_id": "Primary key of the aromatic ring atom.",
        "aromatic_ring_id": "Primary key of the aromatic ring this atom is part of.",
        "atom_id": "Primary key of the atom."
    }
}

comment_on_table_elements(aromatic_ring_atoms, aromatic_ring_atom_comments)

# interactions between aromatic rings
ring_interactions = Table('ring_interactions', metadata,
                          Column('ring_interaction_id', Integer, primary_key=True),
                          Column('biomolecule_id', Integer, nullable=False),
                          Column('aromatic_ring_bgn_id', Integer, nullable=False),
                          Column('aromatic_ring_end_id', Integer, nullable=False),
                          Column('closest_atom_bgn_id', Integer),
                          Column('closest_atom_end_id', Integer),
                          Column('distance', Float(3,2), nullable=False),
                          Column('closest_atom_distance', Float(3,2)),
                          Column('dihedral', Float(5,2), nullable=False),
                          Column('theta', Float(5,2), nullable=False),
                          Column('iota', Float(5,2), nullable=False),
                          Column('interaction_type', String(2)), # "'FF'", "'OF'", "'EE'", "'FT'", "'OT'", "'ET'", "'FE'", "'OE'", "'EF'")
                          schema=schema)

Index('idx_ring_interactions_biomolecule_id', ring_interactions.c.biomolecule_id)
Index('idx_ring_interactions_aromatic_ring_bgn_id', ring_interactions.c.aromatic_ring_bgn_id)
Index('idx_ring_interactions_aromatic_ring_end_id', ring_interactions.c.aromatic_ring_end_id)
Index('idx_ring_interactions_closest_atom_bgn_id', ring_interactions.c.closest_atom_bgn_id)
Index('idx_ring_interactions_closest_atom_end_id', ring_interactions.c.closest_atom_end_id)

ring_interaction_comments = {
    "table": "Contains the parameters of the interaction between two aromatic ring systems.",
    "columns":
    {
        "ring_interaction_id": "Primary key of the ring interaction.",
        "biomolecule_id": "Foreign key of the parent biomolecule.",
        "aromatic_ring_bgn_id": "Primary key of the first aromatic ring.",
        "aromatic_ring_end_id": "Primary key of the second aromatic ring.",
        "closest_atom_bgn_id": "Primary key of the closest atom of the first aromatic ring.",
        "closest_atom_end_id": "Primary key of the closest atom of the second aromatic ring.",
        "distance": "Distance between the aromatic ring centroids.",
        "closest_atom_distance": "Distance between the closest atoms of the aromatic rings.",
        "dihedral": "Dihedral angle between the two planes (normals).",
        "theta": "Signed angle in degrees between the normal of the first aromatic ring and the vector between the two centroids.",
        "iota": "Signed angle in degrees between the normal of the second aromatic ring and the vector between the two centroids.",
        "interaction_type": "Classification of the ring interaction geometry."
    }
}

comment_on_table_elements(ring_interactions, ring_interaction_comments)

atom_ring_interactions = Table('atom_ring_interactions', metadata,
                               Column('atom_ring_interaction_id', Integer, primary_key=True),
                               Column('biomolecule_id', Integer, nullable=False),
                               Column('aromatic_ring_id', Integer, nullable=False),
                               Column('atom_id', Integer, nullable=False),
                               Column('distance', Float(3,2), nullable=False),
                               Column('theta', Float(5,2), nullable=False),
                               Column('interaction_type', String(10)), # CARBONPI, CATIONPI, DONORPI / NULL IF UNKNOWN
                               UniqueConstraint('aromatic_ring_id', 'atom_id', name='unique_interaction'),
                               schema=schema)

Index('idx_atom_ring_interactions_biomolecule_id', atom_ring_interactions.c.biomolecule_id)
Index('idx_atom_ring_interactions_aromatic_ring_id', atom_ring_interactions.c.aromatic_ring_id)
Index('idx_atom_ring_interactions_atom_id', atom_ring_interactions.c.atom_id)

atom_ring_interaction_comments = {
    "table": "Contains the parameters of the interaction between an atom and an aromatic ring system.",
    "columns":
    {
        "atom_ring_interaction_id": "Primary key of the atom-aromatic ring interaction",
        "aromatic_ring_id": "Primary key of the aromating ring participating in this interaction.",
        "atom_id": "Primary key of the atom that is part of this interaction.",
        "distance": "Distance between the atom and the centroid of the aromatic ring.",
        "theta": "Signed angle in degrees between the normal of the aromatic ring and the vector between the centroid and the atom.",
        "interaction_type": "Classification of the atom-aromatic ring interaction."
    }
}

comment_on_table_elements(atom_ring_interactions, atom_ring_interaction_comments)