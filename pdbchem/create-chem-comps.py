import os
import json
import warnings

from openeye.oeiupac import OECreateIUPACName
from openeye.oemolprop import OEGet2dPSA, OEGetXLogP, OEGetFractionCsp3

from eyesopen.oechem import *
from eyesopen.oemolprop import *
from eyesopen.predicates import *
from sqlalchemyutils import OEFP, RDMol, RDBinaryFingerprint

from sqlalchemy import *
from sqlalchemy.event import listen
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects.postgresql import ARRAY

engine      = create_engine(URL(drivername='postgresql+psycopg2', username='adrian',
                                password='1r1d1Um', host='bahamut.bioc.cam.ac.uk',
                                port=5432, database='cryst'),
                            execution_options={'autocommit':True}, echo=False)
connection  = engine.connect()
metadata    = MetaData(bind=connection)

PDBCHEM_SDF_DIR = '/tlbnas/mirror/pdbe/pdbechem/files/sdf'
PDBCHEM_PDB_DIR = '/tlbnas/mirror/pdbe/pdbechem/files/pdb'
SCHEMA          = 'pdbchem'

# CHECK IF RDKIT EXISTS
HAS_RDKIT = connection.execute("SELECT EXISTS (SELECT * FROM pg_catalog.pg_namespace WHERE nspname = 'rdkit')").scalar()

def create_tables():
    '''
    '''
    chem_comps = Table('chem_comps', metadata,
                       Column('chem_comp_id', Integer, primary_key=True),
                       Column('het_id', String(3), nullable=False, index=True),
                       Column('three_letter_code', String(3), index=True),
                       Column('replaced_by_het_id', String(8)),
                       Column('nstd_parent_het_id', Text),
                       Column('subcomponents', ARRAY(String(3))),
                       Column('subcomponent_of_het_ids', ARRAY(String(3))),
                       Column('iupac_name', Text),
                       Column('initial_date', Date),
                       Column('modified_date', Date),
                       Column('ism', Text),
                       Column('inchi', Text),
                       Column('inchikey', Text),
                       Column('mw', Float(8,2)),
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
                       Column('alerts', Integer),
                       Column('qed', Float(4,3)),
                       Column('unwanted_groups', ARRAY(Text)),
                       Column('is_solvent', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_nat_product', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_amino_acid', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_nucleotide', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_saccharide', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_block_buster', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_fragment', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_macrocyclic', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_drug_like', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_lead', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_drug', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_approved_drug', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('is_lig_in_credo', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('has_std_atoms', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('has_butyl', Boolean(create_constraint=False), DefaultClause('false')),
                       Column('has_pains_viol', Boolean(create_constraint=False), DefaultClause('false')),
                       schema=SCHEMA)

    chem_comp_structures = Table('chem_comp_structures', metadata,
                                 Column('het_id', String(3), nullable=False, primary_key=True),
                                 Column('pdb', Text),
                                 Column('sdf', Text),
                                 Column('oeb', LargeBinary),
                                 schema=SCHEMA)

    # table to hold the OpenEye fingerprints
    chem_comp_oefps = Table('chem_comp_oefps', metadata,
                            Column('het_id', String(3), nullable=False, primary_key=True),
                            Column('circular_fp', OEFP),
                            Column('maccs166_fp', OEFP),
                            Column('path_fp', OEFP),
                            Column('tree_fp', OEFP),
                            schema=SCHEMA)

    if HAS_RDKIT:

        chem_comp_rdmols = Table('chem_comp_rdmols', metadata,
                                Column('het_id', String(3), nullable=False, primary_key=True),
                                Column('rdmol', RDMol),
                                schema=SCHEMA)

        chem_comp_rdfps = Table('chem_comp_rdfps', metadata,
                            Column('het_id', String(3), nullable=False, primary_key=True),
                            Column('atompair_fp', RDBinaryFingerprint),
                            Column('torsion_fp', RDBinaryFingerprint),
                            Column('circular_fp', RDBinaryFingerprint),
                            schema=SCHEMA)

    listen(metadata, "after_create", DDL("GRANT SELECT ON ALL TABLES IN SCHEMA pdbchem TO public"))

    metadata.drop_all(checkfirst=True)
    metadata.create_all(checkfirst=True)

    return metadata

def main():
    '''
    '''
    metadata = create_tables()

    chem_comps = metadata.tables['pdbchem.chem_comps']
    chemstruct = metadata.tables['pdbchem.chem_comp_structures']
    chem_comp_oefps = metadata.tables['pdbchem.chem_comp_oefps']
    chem_comp_rdmols = metadata.tables['pdbchem.chem_comp_rdmols']
    chem_comp_rdfps = metadata.tables['pdbchem.chem_comp_rdfps']

    ins_chem_comps = chem_comps.insert()
    ins_chemstruct = chemstruct.insert()

    OEThrow.SetLevel(OEErrorLevel_Error)

    ss_butyl = OESubSearch('[CH2,CH3][CH2][CH2][CH2]')

    ifs = oemolistream()

    dots = OEDots(100,1,'Molecules')

    for filename in os.listdir(PDBCHEM_SDF_DIR):
        name, ext = os.path.splitext(filename)

        if ext != '.sdf': continue

        dots.Update()

        path_sdf = os.path.join(PDBCHEM_SDF_DIR, filename)
        path_pdb = os.path.join(PDBCHEM_PDB_DIR, '{0}.pdb'.format(name))

        ifs.open(str(path_sdf))

        try: chemcomp = ifs.GetOEGraphMols().next()
        except StopIteration: continue

        het_id = chemcomp.GetTitle().split('.')[0]

        OESuppressHydrogens(chemcomp)
        chemcomp = reset_charges(chemcomp)

        # CREATE SENSIBLE IUPAC NAME
        iupac_name = OECreateIUPACName(chemcomp)
        iupac_name = iupac_name if not 'BLAH' in iupac_name else None

        # SMILES
        ism = mol_to_smiles(chemcomp, from3d=True, isomeric=True, reset_charges=True)

        num_hvy_atoms = OECount(chemcomp, OEIsHeavy())
        num_carbons = OECount(chemcomp, OEIsCarbon())
        num_heteroatoms=num_hvy_atoms - num_carbons

        num_ring_systems, ring_atoms = OEDetermineRingSystems(chemcomp)

        try:
            max_ring_size = max([ring_atoms.count(i+1) for i in range(num_ring_systems)])
        except ValueError:
            max_ring_size = 0

        # XLOGP
        try: xlogp = OEGetXLogP(chemcomp)
        except RuntimeError: xlogp = None

        try: qedscore = calc_qed(chemcomp)
        except RuntimeError: qedscore = None

        # CHEMICAL STRUCTURES
        if os.path.exists(path_pdb):
            pdb = open(path_pdb).read()
            pdbmol = mol_from_file(path_pdb)

            try:
                oeb_pdb = mol_to_oeb(pdbmol)
            except ValueError:
                print het_id, path_pdb

        else:
            pdb = None

        sdf = open(path_sdf).read()

        connection.execute(ins_chem_comps,
                           het_id = het_id,
                           three_letter_code = chemcomp.GetTitle()[:3],
                           iupac_name = iupac_name,
                           ism = ism,
                           mw = OECalculateMolecularWeight(chemcomp),
                           num_hvy_atoms=num_hvy_atoms,
                           num_carbons=num_carbons,
                           num_heteroatoms=num_heteroatoms,
                           num_halides = OECount(chemcomp, OEIsHalogen()),
                           het_carb_ratio = num_heteroatoms / float(num_carbons) if num_carbons > 0 else None,
                           formal_charge_count = sum([abs(atom.GetFormalCharge()) for atom in chemcomp.GetAtoms()]),
                           formal_charge_sum = OENetCharge(chemcomp),
                           num_chiral_centers = OECount(chemcomp, OEIsChiralBond()),
                           num_bonds = chemcomp.NumBonds(),
                           num_rotors = OECount(chemcomp, OEIsRotor()),
                           num_ring_systems = num_ring_systems,
                           num_aro_ring_systems = OEDetermineAromaticRingSystems(chemcomp)[0],
                           max_ring_size=max_ring_size,
                           num_lipinski_hbond_acceptors = OECount(chemcomp, OEOrAtom(OEIsNitrogen(), OEIsOxygen())),
                           num_lipinski_hbond_donors = OECount(chemcomp, OEMatchAtom('[#7,#8;!H0]')),
                           tpsa = OEGet2dPSA(chemcomp),
                           xlogp = xlogp,
                           fraction_car = OECount(chemcomp, OEAndAtom(OEIsCarbon(), OEIsAromaticAtom())) / float(num_carbons) if num_carbons > 0 else None,
                           fraction_csp3 = OEGetFractionCsp3(chemcomp),
                           alerts = num_struct_alerts(chemcomp),
                           qed = qedscore,
                           has_std_atoms = 'false' if OECount(chemcomp, OEIsStdAtom()) != chemcomp.NumAtoms() else 'true',
                           has_butyl = 'true' if OECount(chemcomp, ss_butyl) else 'false')

        # INSERT STRUCTURES
        connection.execute(ins_chemstruct,
                           het_id = het_id,
                           pdb = pdb,
                           sdf = sdf,
                           oeb = oeb_pdb)

    dots.Total()

    # substring index on smiles
    DDL("CREATE INDEX idx_chem_comps_ism ON %(fullname)s (SUBSTRING(ism for 64)) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comps)

    # add trigram index on smiles column
    DDL("CREATE INDEX  ON %(fullname)s USING GIST(ism gist_trgm_ops) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comps)

    connection.execute("""
                   INSERT INTO pdbchem.chem_comp_oefps
                   SELECT het_id,
                          openeye.make_circular_fp(ism) AS circular_fp,
                          openeye.make_maccs166_fp(ism) AS maccs166_fp,
                          openeye.make_path_fp(ism) AS path_fp,
                          openeye.make_tree_fp(ism) AS tree_fp
                     FROM pdbchem.chem_comps pc
                   """)

    DDL("CREATE INDEX idx_chem_comp_oefps_circular_fp ON %(fullname)s USING GIST(circular_fp) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comp_oefps)
    DDL("CREATE INDEX idx_chem_comp_oefps_maccs166_fp ON %(fullname)s USING GIST(maccs166_fp) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comp_oefps)
    DDL("CREATE INDEX idx_chem_comp_oefps_path_fp ON %(fullname)s USING GIST(path_fp) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comp_oefps)
    DDL("CREATE INDEX idx_chem_comp_oefps_tree_fp ON %(fullname)s USING GIST(tree_fp) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comp_oefps)

    if HAS_RDKIT:

        # RDKIT CARTRIDGE MOLECULES
        connection.execute("""
                        INSERT   INTO pdbchem.chem_comp_rdmols
                        SELECT   het_id, ism::rdkit.mol
                        FROM     pdbchem.chem_comps
                        WHERE    rdkit.is_valid_smiles(ism::cstring)
                        """)

        DDL("CREATE INDEX idx_chem_comp_rdmols_rdmol ON %(fullname)s USING GIST(rdmol) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comp_rdmols)
        DDL("ALTER TABLE %(fullname)s CLUSTER ON idx_chem_comp_rdmols_rdmol").execute(bind=engine, target=chem_comp_rdmols)

        # CLUSTER TABLE USING THE GIST INDEX
        connection.execute("CLUSTER  pdbchem.chem_comp_rdmols;")

        # RDKIT FINGERPRINTS
        connection.execute("""
                        INSERT   INTO pdbchem.chem_comp_rdfps
                        SELECT   het_id,
                                    rdkit.atompairbv_fp(rdmol) as atompair_fp,
                                    rdkit.torsionbv_fp(rdmol) as torsion_fp,
                                    rdkit.morganbv_fp(rdmol,2) as circular_fp
                        FROM     pdbchem.chem_comp_rdmols
                        """)

        DDL("CREATE INDEX idx_chem_comp_rdfps_atompair_fp ON %(fullname)s USING GIST(atompair_fp) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comp_rdfps)
        DDL("CREATE INDEX idx_chem_comp_rdfps_torsion_fp ON %(fullname)s USING GIST(torsion_fp) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comp_rdfps)
        DDL("CREATE INDEX idx_chem_comp_rdfps_circular_fp ON %(fullname)s USING GIST(circular_fp) WITH (FILLFACTOR=100)").execute(bind=engine, target=chem_comp_rdfps)

    # REPLACED COMPONENTS
    connection.execute("""
                       UPDATE   pdbchem.chem_comps cc
                       SET      inchi = pc.pdbx_inchi,
                                inchikey = pc.pdbx_inchikey,
                                replaced_by_het_id = pc.pdbx_replaced_by,
                                nstd_parent_het_id = pc.mon_nstd_parent_comp_id,
                                subcomponents = pc.pdbx_subcomponent_list
                       FROM     pdbchem.pdb_chem_comps pc
                       WHERE    pc.id = cc.het_id
                       """)

    # AMINO ACIDS
    connection.execute("""
                       UPDATE  pdbchem.chem_comps cc
                       SET     is_amino_acid = true
                       FROM    bio.standard_residues r
                       WHERE   r.three_letter_Code = cc.nstd_parent_het_id
                       """)

    # NUCLEOTIDES
    connection.execute("""
                       UPDATE pdbchem.chem_comps
                          SET is_nucleotide = true
                        WHERE het_id IN ('A','DG','DA','DC','DU','G','C','U','DT','DI','G,A','T','N')
                              OR nstd_parent_het_id IN ('A','DG','DA','DC','DU','G','C','U','DT','DI','G,A','T','N')
                       """)

    # SACCHARIDES
    connection.execute("""
                       UPDATE  pdbchem.chem_comps cc
                       SET     is_saccharide = true
                       FROM    pdbchem.pdb_chem_comps pc
                       WHERE   pc.id = cc.het_id AND UPPER(pc.monomer_type) LIKE '%%SACCHARIDE%%'
                       """)

    # DRUG-LIKE
    connection.execute("""
                       UPDATE  pdbchem.chem_comps
                       SET     is_drug_like = true
                       WHERE   -- REMOVE ALL ENDOGENOUS MOLECULES
                               is_amino_acid = false
                               AND is_saccharide = false
                               AND is_nucleotide = false
                               AND has_std_atoms = true
                               AND has_butyl = false
                               -- LIPINSKI
                               AND num_rotors <= 10
                               AND num_lipinski_hbond_acceptors <= 10
                               AND num_lipinski_hbond_donors <= 5
                               -- GHOSE PROPERTIES
                               AND xlogp BETWEEN -0.4 AND 5.6
                               AND tpsa <= 140
                               AND num_hvy_atoms BETWEEN 20 AND 70
                               AND mw BETWEEN 160 AND 600 -- CHANGED FOR OUTLIERS LIKE IMATINIB
                       """)

    # DRUGS
    connection.execute("""
                       UPDATE   pdbchem.chem_comps cc
                       SET      is_drug = true
                       FROM     pdbchem.chem_comp_drugs cd
                       WHERE    cc.het_id = cd.het_id
                       """)

    # APPROVED DRUGS
    connection.execute("""
                       UPDATE   pdbchem.chem_comps cc
                       SET      is_approved_drug = true
                       FROM     pdbchem.chem_comp_drugs cd
                       WHERE    cc.het_id = cd.het_id AND cd.is_approved = true
                       """)

    #
    connection.execute("""
                       UPDATE   pdbchem.chem_comps cc
                       SET      initial_date = pc.pdbx_initial_date
                       FROM     pdbchem.pdb_chem_comps pc
                       WHERE    cc.het_id = pc.id
                       """)

    connection.execute("""
                       UPDATE   pdbchem.chem_comps cc
                       SET      is_block_buster = openeye.is_block_buster(ism)
                       """)

    connection.execute("""
                       UPDATE   pdbchem.chem_comps cc
                       SET      is_lead = openeye.is_lead(ism)
                       """)

    connection.execute("""
                       UPDATE   pdbchem.chem_comps cc
                       SET      modified_date = pc.pdbx_modified_date
                       FROM     pdbchem.pdb_chem_comps pc
                       WHERE    cc.het_id = pc.id
                       """)

    connection.execute("""
                       UPDATE      pdbchem.chem_comps cc
                        SET         is_solvent = true
                        FROM        pdbchem.astex_solvents a
                        WHERE       a.het_id = cc.het_id;
                       """)

    connection.execute("""
                     UPDATE pdbchem.chem_comps cc
                        SET subcomponent_of_het_ids = sq.het_ids
                       FROM (
                                 WITH subcomponents AS
                                      (
                                       SELECT het_id, unnest(subcomponents) AS subcomponent
                                         FROM pdbchem.chem_comps
                                        WHERE subcomponents IS NOT NULL
                                      )
                               SELECT subcomponent, array_agg(het_id) AS het_ids
                                 FROM subcomponents
                             GROUP BY subcomponent
                            ) sq
                      WHERE sq.subcomponent = cc.het_id
                       """)

    connection.execute("""
                     UPDATE pdbchem.chem_comps cc
                        SET is_lig_in_credo = true
                       FROM credo.ligands l
                      WHERE l.ligand_name = cc.het_id
                       """)

    connection.execute("""
                       UPDATE pdbchem.chem_comps cc
                          SET is_macrocyclic = true
                        WHERE openeye.matches(ism, '[r;!r3;!r4;!r5;!r6;!r7;!r8]')
                       """)

    connection.execute("""
                         UPDATE pdbchem.chem_comps cc
                            SET is_nat_product = true
                           FROM (
                                 WITH intenz as
                                      (
                                       SELECT het_id
                                         FROM (
                                               select unnest(reactants) as chebi_id
                                                 from intenz.reactions
                                                union
                                               select unnest(products) as chebi_id
                                                 from intenz.reactions
                                                union
                                               select unnest(cofactors)
                                                 from intenz.enzymes
                                              ) sq
                                         JOIN scifdw.unichem_pdb_to_chebi u ON u.chebi_id = sq.chebi_id
                                      ),
                                      kegg as
                                      (
                                       SELECT het_id
                                         FROM (
                                               SELECT compound_id FROM kegg.enzyme_products
                                                UNION
                                               SELECT compound_id FROM kegg.enzyme_substrates
                                              ) sq
                                         JOIN scifdw.unichem_pdb_to_kegg u ON u.compound_id = sq.compound_id
                                      )
                                 SELECT * FROM intenz UNION SELECT * FROM kegg
                                ) sq
                          WHERE sq.het_id = cc.het_id
                       """)


    connection.execute("""
                         UPDATE pdbchem.chem_comps cc
                            SET unwanted_groups = sq.groups
                           FROM (
                                   SELECT het_id, array_agg(DISTINCT group_name) AS groups
                                     FROM pdbchem.chem_comps cc, chemistry.unwanted_groups u
                                    WHERE openeye.matches(ism, u.smarts)
                                 GROUP BY het_id
                                ) sq
                          WHERE sq.het_id = cc.het_id;
                       """)

    connection.execute("""
                        DO $$
                            DECLARE
                                chem_id INTEGER;
                            BEGIN
                                FOR chem_id IN SELECT chem_comp_id FROM pdbchem.chem_comps cc ORDER BY 1
                                LOOP
                                    EXECUTE
                                    '
                                     UPDATE pdbchem.chem_comps cc
                                        SET has_pains_viol = true
                                       FROM chemistry.pains p
                                      WHERE cc.chem_comp_id = $1
                                            AND openeye.matches(ism, smarts);
                                    ' USING chem_id;

                                    RAISE NOTICE 'Updated PAINS for chemical component %%', chem_id;
                                END LOOP;
                        END$$;
                       """)

    # use a basic rule of three
    connection.execute("""
                       UPDATE pdbchem.chem_comps
                          SET is_fragment = true
                        WHERE mw between 80 and 300
                              AND has_std_atoms = true
                            AND num_carbons > 0
                            AND num_heteroatoms > 0
                            AND num_halides <= 4
                            AND het_carb_ratio < 1.75
                            AND num_lipinski_hbond_acceptors BETWEEN 1 AND 6
                            AND num_lipinski_hbond_donors <= 3
                            AND num_rotors <= 3
                            AND num_hvy_atoms BETWEEN 6 AND 18
                            AND num_ring_systems BETWEEN 1 AND 3
                            AND max_ring_size <= 10
                            -- total number of aromatic atoms:
                            -- maximum 2 naphtalene + disconnected benzene
                            AND openeye.count_match_atoms(ism, '[a]') <= 15
                            AND openeye.count_match_atoms(ism, '[#35,#53]') < 2
                            AND xlogp BETWEEN -4 AND 4
                            AND NOT (max_ring_size = num_hvy_atoms AND num_hvy_atoms <= 6)
                            AND has_butyl = false
                            -- no PAINS compounds
                            AND has_pains_viol = False
                            -- no saccharides
                            AND is_saccharide = false
                            AND NOT (openeye.matches(ism, 'O1[CR2][CR2][CR2][CR2]1')
                                     OR openeye.matches(ism,'O1[CR2][CR2][CR2][CR2][CR2]1'))
                       """)

main()
