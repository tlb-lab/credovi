#!/usr/bin/python
import os
import json

from openeye.oechem import *
from openeye.oequacpac import OESetNeutralpHModel
from openeye.oeiupac import OECreateIUPACName
from openeye.oemolprop import OEGet2dPSA, OEGetXLogP, OEGetFractionCsp3
from sqlalchemy import *
from sqlalchemy.event import listen
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects.postgresql import ARRAY

from eyesopen.oechem import mol_to_smiles, reset_charges
from eyesopen.predicates import OEIsStdAtom
from eyesopen.oerecap import decompose

engine      = create_engine(URL(drivername='postgresql+psycopg2', username='adrian',
                                password='1r1d1Um', host='bahamut.bioc.cam.ac.uk',
                                port=5432, database='cryst'),
                            execution_options={'autocommit':True}, echo=False)
connection  = engine.connect()
metadata    = MetaData(bind=connection)

PDBCHEM_PDB_DIR = '/tlbnas/mirror/pdbe/pdbechem/files/pdb'
PDBCHEM_SDF_DIR = '/tlbnas/mirror/pdbe/pdbechem/files/sdf'
RAW_FRAG_FILE   = '/tlbnas/temp/bahamut/raw_fragments.pdbchem'

def create_tables():
    """
    """
    raw_fragments = Table('raw_fragments', metadata,
                          Column('het_id', String(15), index=True),
                          Column('order_parent', Integer()),
                          Column('ism_parent', Text),
                          Column('parent_atom_names', Text(), nullable=False),
                          Column('order_child', Integer()),
                          Column('ism_child', Text),
                          Column('child_atom_names', Text()),
                          schema='pdbchem', prefixes = ['unlogged'])

    # CUSTOM INDEX ON SMILES
    listen(raw_fragments, "after_create",
           DDL("CREATE INDEX idx_raw_fragments_ism_parent ON %(fullname)s (substring(ism_parent for 64));"))
    listen(raw_fragments, "after_create",
           DDL("CREATE INDEX idx_raw_fragments_ism_child ON %(fullname)s (substring(ism_child for 64));"))

    fragments = Table('fragments', metadata,
                      Column('fragment_id', Integer, primary_key=True),
                      Column('ism', Text, nullable=False),
                      Column('iupac_name', Text),
                      Column('mw', Float(8,2)),
                      Column('num_atoms', Integer),
                      Column('num_hvy_atoms', Integer),
                      Column('num_carbons', Integer),
                      Column('num_heteroatoms', Integer),
                      Column('num_halides', Integer),
                      Column('het_carb_ratio', Float(4,2)),
                      Column('formal_charge_count', Integer),
                      Column('formal_charge_sum', Integer),
                      Column('num_chiral_centers', Integer),
                      Column('num_bonds', Integer),
                      Column('num_rotors', Integer),
                      Column('num_ring_systems', Integer),
                      Column('num_aro_ring_systems', Integer),
                      Column('max_ring_size', Integer),
                      Column('num_lipinski_hbond_acceptors', Integer),
                      Column('num_lipinski_hbond_donors', Integer),
                      Column('tpsa', Float(6,2)),
                      Column('xlogp', Float(4,2)),
                      Column('fraction_car', Float(4,3)),
                      Column('fraction_csp3', Float(4,3)),
                      Column('unwanted_groups', ARRAY(Text)),
                      Column('is_terminal', Boolean(create_constraint=False), DefaultClause('false')),
                      Column('is_shared', Boolean(create_constraint=False), DefaultClause('false')),
                      Column('is_unfragmented', Boolean(create_constraint=False), DefaultClause('false')),
                      Column('is_rule_of_three', Boolean(create_constraint=False), DefaultClause('false')),
                      Column('is_purchasable', Boolean(create_constraint=False), DefaultClause('false')),
                      Column('is_true_pdb_fragment', Boolean(create_constraint=False), DefaultClause('false')),
                      Column('has_std_atoms', Boolean(create_constraint=False), DefaultClause('false')),
                      Column('has_butyl', Boolean(create_constraint=False), DefaultClause('false')),
                      Column('has_pains_viol', Boolean(create_constraint=False), DefaultClause('false')),
                      schema='pdbchem'
                      )

    # CUSTOM INDEX ON SMILES
    listen(fragments, "after_create", DDL("CREATE INDEX idx_fragments_ism ON %(fullname)s (substring(ism for 64));"))

    hierarchies = Table('fragment_hierarchies', metadata,
                        Column('fragment_hierarchy_id', Integer, primary_key=True),
                        Column('het_id', String(15), nullable=False, index=True),
                        Column('parent_id', Integer(), ForeignKey(fragments.c.fragment_id), nullable=False),
                        Column('order_parent', Integer),
                        Column('child_id', Integer(), ForeignKey(fragments.c.fragment_id)),
                        Column('order_child', Integer),
                        Column('is_smallest_parent', Boolean(create_constraint=False), DefaultClause('false')),
                        schema='pdbchem'
    )

    Index('idx_fragment_hierarchies_parent', hierarchies.c.parent_id, hierarchies.c.order_parent)
    Index('idx_fragment_hierarchies_child', hierarchies.c.child_id, hierarchies.c.order_child)

    chem_comp_fragments = Table('chem_comp_fragments', metadata,
                                Column('chem_comp_fragment_id', Integer, primary_key=True),
                                Column('het_id', String(15), nullable=False),
                                Column('fragment_id', Integer(), ForeignKey(fragments.c.fragment_id), nullable=False, index=True),
                                Column('is_root', Boolean(create_constraint=False), DefaultClause('false')),
                                schema='pdbchem')

    Index('idx_chem_comp_fragments_het_id', chem_comp_fragments.c.het_id, chem_comp_fragments.c.fragment_id, unique=True)

    chem_comp_fragment_atoms = Table('chem_comp_fragment_atoms', metadata,
                                     Column('chem_comp_fragment_atom_id', Integer, primary_key=True),
                                     Column('chem_comp_fragment_id', Integer, ForeignKey(chem_comp_fragments.c.chem_comp_fragment_id), nullable=False, index=True),
                                     Column('hit', Integer,nullable=False),
                                     Column('pdb_name', String(4)),
                                     schema='pdbchem')

    metadata.drop_all(checkfirst=True)
    metadata.create_all(checkfirst=True)

    return metadata

def clean_node_ism(fragment):
    """
    """
    return mol_to_smiles(fragment, isomeric=True, from3d=True, reset_charges=True).replace('\\','\\\\')

def traverse(node, order, het_id):
    """
    Traverse the hierarchies dictionary in recursive fashion to get all parent->child relations
    """
    raw = open(RAW_FRAG_FILE,'a')

    # TO KEEP COMPONENTS IN THE HIERARCHY THAT CANNOT BE FRAGMENTED (ADENINE)
    if (order == 0) and not len(node.children.keys()):
        ism_parent = clean_node_ism(node.mol)
        string = '\t'.join([het_id, str(order), ism_parent, node.atoms, '\N', '\N', '\N'])

        print >> raw, string
        raw.flush()

    else:
        while len(node.children.keys()):
            ism_parent = clean_node_ism(node.mol)
            ism_child, child_node = node.children.popitem()
            ism_child = clean_node_ism(child_node.mol)

            string='\t'.join([het_id, str(order), ism_parent, node.atoms, str(order + 1), ism_child, child_node.atoms])

            print >> raw, string

            raw.flush()

            traverse(child_node, order+1, het_id)

    raw.close()

def main():
    """
    """
    elements = (OEElemNo_He,OEElemNo_Hf,OEElemNo_Hg,OEElemNo_Ho,OEElemNo_Hs)

    pdbifs = oemolistream()
    pdbifs.SetFormat(OEFormat_PDB)
    sdfifs = oemolistream()
    sdfifs.SetFormat(OEFormat_SDF)

    # SUPPRESS WARNINGS OF ZERO REACTANTS
    OEThrow.SetLevel(OEErrorLevel_Error)

    #
    if os.path.exists(RAW_FRAG_FILE): os.unlink(RAW_FRAG_FILE)

    dots = OEDots(100, 1, "Chemical components")

    for chem_pdb_file in os.listdir(PDBCHEM_PDB_DIR):
        dots.Update()

        het_id = chem_pdb_file[:-4]

        pdb_path = os.path.join(PDBCHEM_PDB_DIR, chem_pdb_file)
        sdf_path = os.path.join(PDBCHEM_SDF_DIR, '{}.sdf'.format(het_id))

        if not os.path.exists(pdb_path) or not os.path.exists(sdf_path):
            continue

        pdbifs.open(str(pdb_path))
        sdfifs.open(str(sdf_path))

        chemcomp, pdb = OEGraphMol(), OEGraphMol()

        OEReadMolecule(sdfifs, chemcomp)
        OEReadMolecule(pdbifs, pdb)

        OESuppressHydrogens(chemcomp)
        OESuppressHydrogens(pdb)

        # copy the PDB atom names to the SDF molecule
        for i,j in zip(chemcomp.GetAtoms(), pdb.GetAtoms()):
            i.SetName(j.GetName())

        # IGNORE STRUCTURES WITH WRONGLY PERCEIVED HYDROGENS
        if any(atom.GetAtomicNum() in elements for atom in chemcomp.GetAtoms()):
            continue

        # IGNORE VERY LARGE MOLECULES AND SMALL IONS/FRAGMENTS
        if chemcomp.NumAtoms() > 100 or chemcomp.NumAtoms() < 3: continue

        # RESET FORMAL CHARGES
        chemcomp = reset_charges(chemcomp)

        node = decompose(chemcomp, min_fragment_size=0, allow_butyl=True)

        # TRAVERSE NODE AND WRITE DETAILS TO FILE
        traverse(node, 0, het_id)

    dots.Total()

    metadata = create_tables()

    DDL("COPY pdbchem.raw_fragments FROM '{}'".format(RAW_FRAG_FILE)).execute(bind=engine)

    connection.execute("GRANT SELECT ON ALL TABLES IN SCHEMA pdbchem TO PUBLIC")

    connection.execute("""
                       DELETE FROM pdbchem.raw_fragments rf
                        USING pdbchem.chem_comps cc
                        WHERE cc.het_id = rf.het_id
                              AND cc.replaced_by_het_id IS NOT NULL
                       """)

    # CREATE UNIQUE REPRESENTATION OF FRAGMENTS ORDERED BY THE NUMBER OF HEAVY ATOMS
    connection.execute("""
                       INSERT       INTO pdbchem.fragments(ism, num_atoms)
                       SELECT       ism, num_atoms
                       FROM         (
                                    SELECT      ism, openeye.smiles_atom_count(ism) as num_atoms
                                    FROM        (
                                                SELECT  openeye.standardise_smiles(ism_parent, true, true) as ism
                                                FROM    pdbchem.raw_fragments rf
                                                ) parents
                                    UNION
                                    SELECT      ism, openeye.smiles_atom_count(ism) as num_atoms
                                    FROM        (
                                                SELECT  openeye.standardise_smiles(ism_child, true, true) as ism
                                                FROM    pdbchem.raw_fragments rf
                                                WHERE   ism_child IS NOT NULL
                                                ) children
                                    ORDER BY    2
                                    ) fragments;
                       """)

    # CREATE THE FRAGMENT HIERARCHY
    connection.execute("""
                       INSERT INTO pdbchem.fragment_hierarchies(het_id, parent_id, order_parent, child_id, order_child)
                       SELECT      DISTINCT r.het_id, a.fragment_id, order_parent, b.fragment_id, order_child
                       FROM        pdbchem.raw_fragments r
                       JOIN        pdbchem.fragments a ON openeye.standardise_smiles(r.ism_parent, true, true) = a.ism
                       LEFT JOIN   pdbchem.fragments b ON openeye.standardise_smiles(r.ism_child, true, true) = b.ism
                       ORDER BY    het_id, order_parent, a.fragment_id;
                       """)


    # LINK CHEMICAL COMPONENTS TO FRAGMENTS
    connection.execute("""
                        INSERT      INTO pdbchem.chem_comp_fragments(het_id, fragment_id)
                        SELECT      DISTINCT het_id, fh.parent_id as fragment_id
                        FROM        pdbchem.fragment_hierarchies fh
                        UNION
                        SELECT      DISTINCT het_id, fh.child_id as fragment_id
                        FROM        pdbchem.fragment_hierarchies fh
                        WHERE       child_id IS NOT NULL
                        ORDER BY    1,2;
                       """)

    connection.execute("""
                        INSERT      INTO pdbchem.chem_comp_fragment_atoms(chem_comp_fragment_id, hit, pdb_name)
                        SELECT      cf.chem_comp_fragment_id,
                                    ROW_NUMBER() OVER (PARTITION BY cf.chem_comp_fragment_id) AS hit,
                                    unnest(string_to_array(sq.atom_names,'+'))
                        FROM        pdbchem.fragments f
                        JOIN        pdbchem.chem_comp_fragments cf
                                    ON cf.fragment_id = f.fragment_id
                        JOIN        (
                                    SELECT  DISTINCT het_id, openeye.standardise_smiles(ism_parent, true, true) AS ism, parent_atom_names AS atom_names
                                    FROM    pdbchem.raw_fragments rf
                                    UNION
                                    SELECT  DISTINCT het_id, openeye.standardise_smiles(ism_child, true, true), child_atom_names
                                    FROM    pdbchem.raw_fragments rf
                                    ) sq
                                    ON sq.ism = f.ism
                                    AND sq.het_id = cf.het_id
                        ORDER BY    1,2
                       """)

    # SET A FLAG FOR TERMINAL FRAGMENTS
    connection.execute("""
                       UPDATE  pdbchem.fragments f
                       SET     is_terminal = true
                       FROM    (
                                SELECT      DISTINCT f.fragment_id
                                FROM        pdbchem.fragments f
                                LEFT JOIN   pdbchem.fragment_hierarchies h ON h.parent_id = f.fragment_id
                                WHERE       h.parent_id IS NULL
                       ) sq
                       WHERE   f.fragment_id = sq.fragment_id
                       """)

    connection.execute("""
                        UPDATE  pdbchem.fragments f
                        SET     is_terminal = true
                        FROM    pdbchem.fragment_hierarchies h
                        WHERE   h.parent_id = f.fragment_id
                                AND h.child_id IS NULL
                       """)

    # FLAG FOR FRAGMENTS THAT CANNOT BE FRAGMENTED
    connection.execute("""
                        UPDATE  pdbchem.fragments f
                        SET     is_unfragmented = true
                        FROM    (
                                SELECT      DISTINCT f.fragment_id
                                FROM        pdbchem.fragments f
                                JOIN        pdbchem.Fragment_hierarchies h1 ON f.fragment_id = h1.parent_id
                                LEFT JOIN   pdbchem.Fragment_hierarchies h2 ON h1.parent_id = h2.child_id
                                WHERE       h1.child_id IS NULL
                                            AND h2.child_id IS NULL
                                ) sq
                        WHERE   sq.fragment_id = f.fragment_id
                       """)

    # FLAG FRAGMENTS THAT ARE SHARED BY MORE THAN ONE CHEMICAL COMPONENT
    connection.execute("""
                        UPDATE  pdbchem.fragments f
                        SET     is_shared = true
                        FROM    (
                                SELECT      fragment_id, COUNT(DISTINCT cf.het_id) het_ids
                                FROM        pdbchem.chem_comp_fragments cf
                                JOIN        pdbchem.chem_comps cc ON cc.het_id = cf.het_id
                                WHERE       cc.replaced_by_het_id IS NULL
                                            AND LENGTH(cc.het_id) <= 3
                                GROUP BY    fragment_id
                                HAVING      COUNT(DISTINCT cf.het_id)  > 1
                                ) sq
                        WHERE   sq.fragment_id = f.fragment_id
                       """)

    connection.execute("""
                        UPDATE pdbchem.chem_comp_fragments cf
                           SET is_root = true
                          FROM (
                                  SELECT het_id, max(fragment_id) as fragment_id
                                    FROM pdbchem.chem_comp_fragments
                                GROUP BY het_id
                               ) sq
                         WHERE cf.het_id = sq.het_id AND cf.fragment_id = sq.fragment_id
                       """)

    connection.execute("""
                         UPDATE pdbchem.fragments f
                            SET is_true_pdb_fragment = true
                           FROM (
                                 SELECT DISTINCT cf.fragment_id
                                   FROM pdbchem.chem_comps cc
                                   JOIN pdbchem.chem_comp_fragments cf USING(het_id)
                                  WHERE cc.is_fragment = true AND cf.is_root = true
                                ) sq
                          WHERE sq.fragment_id = f.fragment_id
                       """)

    connection.execute("""
                         UPDATE pdbchem.fragments f
                            SET is_purchasable = true
                           FROM emolecules.molecules m
                          WHERE openeye.standardise_smiles(f.ism) = m.ism;
                       """)

    connection.execute("""
                         UPDATE pdbchem.fragments f
                            SET is_purchasable = true
                           FROM zinc.purchasable z
                          WHERE openeye.standardise_smiles(f.ism) = z.ism;
                       """)

    connection.execute("""
                        DO $$
                            DECLARE
                                frag_id INTEGER;
                            BEGIN
                                FOR frag_id IN SELECT fragment_id FROM pdbchem.fragments f ORDER BY 1
                                LOOP
                                    EXECUTE
                                    '
                                     UPDATE pdbchem.fragments f
                                        SET has_pains_viol = true
                                       FROM chemistry.pains p
                                      WHERE f.fragment_id = $1
                                            AND openeye.matches(ism, smarts);
                                    ' USING frag_id;

                                    RAISE NOTICE 'Updated PAINS for fragment %%', frag_id;
                                END LOOP;
                        END$$;
                       """)

    connection.execute("""
                        DO $$
                            DECLARE
                                frag_id INTEGER;
                            BEGIN
                                FOR frag_id IN SELECT fragment_id FROM pdbchem.fragments f ORDER BY 1
                                LOOP
                                    EXECUTE
                                    '
                                     UPDATE pdbchem.fragments f
                                        SET unwanted_groups = sq.groups
                                       FROM (
                                               SELECT fragment_id, array_agg(DISTINCT group_name) AS groups
                                                 FROM pdbchem.fragments f, chemistry.unwanted_groups u
                                                WHERE fragment_id = $1
                                                      AND openeye.matches(ism, u.smarts)
                                             GROUP BY fragment_id
                                            ) sq
                                      WHERE sq.fragment_id = f.fragment_id;
                                    ' USING frag_id;

                                    RAISE NOTICE 'Updated unwanted groups for fragment %%', frag_id;
                                END LOOP;
                        END$$;
                       """)

    # update fragment properties

    fragments = Table('fragments', metadata, schema='pdbchem', autoload=True)
    update = fragments.update()

    OEThrow.SetLevel(OEErrorLevel_Error)
    ss_butyl = OESubSearch('[CH2][CH2][CH2][CH2]')

    fragment = OEGraphMol()
    statement = select([fragments.c.fragment_id,fragments.c.ism])

    dots = OEDots(100, 1, "Fragments")

    for fragment_id, ism in connection.execute(statement).fetchall():
        dots.Update()

        OEParseSmiles(fragment, str(ism))

        iupac_name = OECreateIUPACName(fragment)
        iupac_name = iupac_name if not 'BLAH' in iupac_name else None

        # XLOGP
        # does not work with formal charges, has to be neutralized!
        try: xlogp = OEGetXLogP(fragment)
        except RuntimeError: xlogp = None

        num_hvy_atoms = OECount(fragment, OEIsHeavy())
        num_carbons = OECount(fragment, OEIsCarbon())
        num_heteroatoms = num_hvy_atoms - num_carbons
        num_car = OECount(fragment, OEAndAtom(OEIsCarbon(), OEIsAromaticAtom()))

        num_ring_systems, ring_atoms = OEDetermineRingSystems(fragment)

        try:
            max_ring_size = max([ring_atoms.count(i+1) for i in range(num_ring_systems)])
        except ValueError:
            max_ring_size = 0

        values = dict(iupac_name = iupac_name,
                      mw = OECalculateMolecularWeight(fragment),
                      num_hvy_atoms = num_hvy_atoms,
                      num_carbons = num_carbons,
                      num_heteroatoms = num_heteroatoms,
                      num_halides = OECount(fragment, OEIsHalogen()),
                      het_carb_ratio = num_heteroatoms / float(num_carbons) if num_carbons > 0 else None,
                      formal_charge_count = sum([abs(atom.GetFormalCharge()) for atom in fragment.GetAtoms()]),
                      formal_charge_sum = OENetCharge(fragment),
                      num_chiral_centers = OECount(fragment, OEIsChiralBond()),
                      num_bonds = fragment.NumBonds(),
                      num_rotors = OECount(fragment, OEIsRotor()),
                      num_ring_systems=num_ring_systems,
                      num_aro_ring_systems = OEDetermineAromaticRingSystems(fragment)[0],
                      max_ring_size=max_ring_size,
                      num_lipinski_hbond_acceptors = OECount(fragment, OEOrAtom(OEIsNitrogen(), OEIsOxygen())),
                      num_lipinski_hbond_donors = OECount(fragment, OEMatchAtom('[#7,#8;!H0]')),
                      tpsa = OEGet2dPSA(fragment),
                      xlogp = xlogp,
                      fraction_car = num_car / float(num_carbons) if num_carbons > 0 else None,
                      fraction_csp3 = OEGetFractionCsp3(fragment) if num_carbons > 0 else None,
                      has_std_atoms = 'false' if OECount(fragment, OEIsStdAtom()) != fragment.NumAtoms() else 'true',
                      has_butyl = 'true' if OECount(fragment, ss_butyl) else 'false'
                      )

        connection.execute(update.where(fragments.c.fragment_id==fragment_id)
                           .values(**values))

        fragment.Clear()

    dots.Total()

    connection.execute("""
                    UPDATE pdbchem.fragments set is_rule_of_three = true
                     WHERE num_lipinski_hbond_donors <= 3
                       AND num_lipinski_hbond_acceptors <= 3
                       AND xlogp <= 3
                       AND num_rotors <= 3
                       AND tpsa < 60
                       AND mw between 60 AND 300
                       AND num_hvy_atoms >= 5;
                       """)

main()
