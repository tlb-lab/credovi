DROP    TABLE IF EXISTS credo.ligand_fragments;
CREATE  TABLE credo.ligand_fragments (
        ligand_fragment_id SERIAL NOT NULL,
        biomolecule_id INTEGER NOT NULL,
        ligand_id INTEGER NOT NULL,
        ligand_component_id INTEGER NOT NULL,
        fragment_id INTEGER NOT NULL,
        hit INTEGER NOT NULL,
        fcd FLOAT DEFAULT NULL,
        is_root boolean NOT NULL DEFAULT FALSE,
        is_interacting boolean NOT NULL DEFAULT FALSE,
        CONSTRAINT ligand_fragments_pkey PRIMARY KEY (ligand_fragment_id)
);

CREATE  INDEX idx_ligand_fragments_ligand_id
        ON credo.ligand_fragments
        USING btree (ligand_id);

CREATE  INDEX idx_ligand_fragments_ligand_component_id
        ON credo.ligand_fragments
        USING btree (ligand_component_id);

CREATE  INDEX idx_ligand_fragments_fragment_id
        ON credo.ligand_fragments
        USING btree (fragment_id);

GRANT   SELECT ON TABLE credo.ligand_fragments TO credouser;

DROP    TABLE IF EXISTS credo.ligand_fragment_atoms;
CREATE  TABLE credo.ligand_fragment_atoms (
        ligand_fragment_atom_id SERIAL NOT NULL,
        ligand_fragment_id INTEGER NOT NULL,
        atom_id INTEGER NOT NULL,
        CONSTRAINT ligand_fragment_atoms_pkey PRIMARY KEY (ligand_fragment_atom_id)
);

CREATE  INDEX idx_ligand_fragment_atoms_ligand_fragment_id
        ON credo.ligand_fragment_atoms
        USING btree (ligand_fragment_id);

CREATE INDEX idx_ligand_fragment_atoms_atom_id
        ON credo.ligand_fragment_atoms
        USING btree (atom_id);

GRANT   SELECT ON TABLE credo.ligand_fragment_atoms TO credouser;

DROP    TABLE IF EXISTS credo.ligand_fcd;
CREATE  TABLE credo.ligand_fcd (
        ligand_fragment_id INTEGER NOT NULL,
        num_int_atoms INTEGER DEFAULT NULL,
        num_covalent INTEGER NOT NULL,
        num_vdw_clash INTEGER NOT NULL,
        num_vdw INTEGER NOT NULL,
        num_proximal INTEGER NOT NULL,
        num_hbond INTEGER NOT NULL,
        num_weak_hbond INTEGER NOT NULL,
        num_xbond INTEGER NOT NULL,
        num_ionic INTEGER NOT NULL,
        num_metal_complex INTEGER NOT NULL,
        num_aromatic INTEGER NOT NULL,
        num_hydrophobic INTEGER NOT NULL,
        num_carbonyl INTEGER NOT NULL,
        fcd FLOAT DEFAULT NULL,
        CONSTRAINT ligand_fcd_pkey PRIMARY KEY (ligand_fragment_id)
);

GRANT   SELECT ON TABLE credo.ligand_fcd TO credopvuser;

-- CREATE A MAPPING BETWEEN LIGANDS AND FRAGMENTS
DO $$
    DECLARE
        lig_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL LIGANDS THAT DO NOT HAVE CLASHES AND ARE NOT DISORDERED
        FOR lig_id IN SELECT ligand_id FROM credo.ligands WHERE is_disordered = false AND is_clashing = false ORDER BY 1
        LOOP
            -- INSERT LIGAND FRAGMENT ATOMS
            EXECUTE
            '
            INSERT      INTO credo.ligand_fragments(biomolecule_id, ligand_id, ligand_component_id, fragment_id, hit)
            SELECT      DISTINCT r.biomolecule_id, lc.ligand_id, lc.ligand_component_id, fragment_id, hit
            FROM        credo.ligand_components lc
            JOIN        credo.residues r ON r.residue_id = lc.residue_id
            JOIN        pdbchem.chem_comp_fragments cf ON r.res_name = cf.het_id
            JOIN        pdbchem.chem_comp_fragment_atoms cfa ON cfa.chem_comp_fragment_id = cf.chem_comp_fragment_id
            WHERE       lc.ligand_id = $1
            ORDER BY    1,2,3;
            ' USING lig_id;

            RAISE NOTICE 'inserted ligand fragments for ligand %', lig_id;
        END LOOP;
END$$;

-- SET A FLAG OF ROOT LIGAND FRAGMENTS
UPDATE credo.ligand_fragments lf
   SET is_root = true
  FROM credo.ligand_components lc,
       pdbchem.fragment_hierarchies fh
 WHERE lc.ligand_component_id = lf.ligand_component_id
       AND fh.parent_id = lf.fragment_id
       AND fh.order_parent = 0
       AND fh.het_id = lc.het_id;

-- CREATE A MAPPING BETWEEN LIGAND FRAGMENTS AND THEIR ATOMS
DO $$
    DECLARE
        lig_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL LIGANDS THAT DO NOT HAVE CLASHES AND ARE NOT DISORDERED
        FOR lig_id, biomol_id IN   SELECT DISTINCT l.ligand_id, l.biomolecule_id
                                     FROM credo.ligands l
                                     JOIN credo.ligand_fragments USING(ligand_id)
                                 ORDER BY 1
        LOOP
            -- INSERT LIGAND FRAGMENT ATOMS
            EXECUTE
            '
               INSERT INTO credo.ligand_fragment_atoms(ligand_fragment_id, atom_id)
               SELECT DISTINCT lf.ligand_fragment_id, a.atom_id
                 FROM credo.ligand_fragments lf
                      -- GET THE LIGAND COMPONENT ATOMS
                 JOIN credo.ligand_components lc ON lc.ligand_component_id = lf.ligand_component_id
                 JOIN credo.residues r ON r.residue_id = lc.residue_id
                 JOIN credo.atoms a on a.residue_id = lc.residue_id
                      -- LINK THEM TO THE FRAGMENT ATOMS
                 JOIN pdbchem.chem_comp_fragments cf
                      ON cf.het_id = r.res_name AND cf.fragment_id = lf.fragment_id
                 JOIN pdbchem.chem_comp_fragment_atoms cfa
                      ON cfa.chem_comp_fragment_id = cf.chem_comp_fragment_id
                      AND lf.hit = cfa.hit
                      AND cfa.pdb_name = a.atom_name
                WHERE lc.ligand_id = $1 AND a.biomolecule_id = $2
             ORDER BY 1,2;
            ' USING lig_id, biomol_id;

            RAISE NOTICE 'inserted ligand fragment atoms for ligand %', lig_id;
        END LOOP;
END$$;

-- FRAGMENT CONTACT DENSITY FOR EACH LIGAND FRAGMENT
DO $$
    DECLARE
        lig_frag_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR lig_frag_id, biomol_id IN   SELECT lf.ligand_fragment_id, lf.biomolecule_id
                                          FROM credo.ligand_fragments lf
                                          JOIN credo.ligands USING(ligand_id)
                                      ORDER BY 1
        LOOP
            EXECUTE
                '
                INSERT      INTO credo.ligand_fcd(ligand_fragment_id, num_int_atoms,
                                                  num_covalent, num_vdw_clash, num_vdw, num_proximal,
                                                  num_hbond, num_weak_hbond, num_xbond, num_ionic,
                                                  num_metal_complex, num_aromatic, num_hydrophobic,
                                                  num_carbonyl)
                WITH        fragment_contacts AS
                            (
                            SELECT      lfa.ligand_fragment_id,
                                        cs.atom_bgn_id as atom_fragment_id, cs.atom_end_id as atom_id, cs.*
                            FROM        credo.ligand_fragment_atoms lfa
                                        -- FRAGMENT ATOM IS BGN
                            JOIN        credo.contacts cs ON cs.atom_bgn_id = lfa.atom_id
                            WHERE       lfa.ligand_fragment_id = $1
                                        AND cs.biomolecule_id = $2
                                        AND cs.is_same_entity = false
                                        AND cs.distance <= 4.5
                                        -- INTERACTION MUST BE WITH POLYMER ATOM (PROT/DNA/RNA)
                                        AND cs.structural_interaction_type_bm & 56 > 0
                            UNION ALL
                            SELECT      lfa.ligand_fragment_id,
                                        cs.atom_end_id as atom_fragment_id, cs.atom_bgn_id as atom_id, cs.*
                            FROM        credo.ligand_fragment_atoms lfa
                                        -- FRAGMENT ATOM IS END
                            JOIN        credo.contacts cs ON cs.atom_end_id = lfa.atom_id
                            WHERE       lfa.ligand_fragment_id = $1
                                        AND cs.biomolecule_id = $2
                                        AND cs.is_same_entity = false
                                        AND cs.distance <= 4.5
                                        -- INTERACTION MUST BE WITH POLYMER ATOM (PROT/DNA/RNA)
                                        AND cs.structural_interaction_type_bm & 3584 > 0
                            )
                SELECT      fc.ligand_fragment_id,
                            COUNT(DISTINCT fc.atom_id) AS num_int_atoms,
                            SUM(fc.is_covalent::int) as num_covalent,
                            SUM(fc.is_vdw_clash::int) as num_vdw_clash,
                            SUM(fc.is_vdw::int) AS num_vdw,
                            SUM(fc.is_proximal::int) AS num_proximal,
                            SUM(fc.is_hbond::int) AS num_hbond,
                            SUM(fc.is_weak_hbond::int) AS num_weak_hbond,
                            SUM(fc.is_xbond::int) as num_xbond,
                            SUM(fc.is_ionic::int) as num_ionic,
                            SUM(fc.is_metal_complex::int) AS num_metal_complex,
                            SUM(fc.is_aromatic::int) AS num_aromatic,
                            SUM(fc.is_hydrophobic::int) AS num_hydrophobic,
                            SUM(fc.is_carbonyl::int) AS num_carbonyl
                FROM        fragment_contacts fc
                GROUP BY    fc.ligand_fragment_id
                ORDER BY    fc.ligand_fragment_id
                ' USING lig_frag_id, biomol_id;

                RAISE NOTICE 'calculated fragment contact densities for ligand fragment %', lig_frag_id;
        END LOOP;
END$$;

-- FRAGMENT CONTACT DENSITY FOR EACH LIGAND FRAGMENT
DO $$
    DECLARE
        lig_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR lig_id, biomol_id IN   SELECT DISTINCT lf.ligand_id, lf.biomolecule_id
                                     FROM credo.ligand_fragments lf
                                     JOIN credo.ligands USING(ligand_id)
                                 ORDER BY 1
        LOOP
            EXECUTE
                '
                UPDATE credo.ligand_fcd f
                   SET fcd = sq.fcd
                  FROM (
                        WITH W AS
                        (
                           SELECT sq.ligand_id, COUNT(DISTINCT sq.atom_id) AS num_int_atoms
                             FROM (
                                   SELECT DISTINCT h.ligand_id, cs.atom_end_id as atom_id
                                     FROM credo.hetatms h
                                     JOIN credo.contacts cs ON cs.atom_bgn_id = h.atom_id
                                    WHERE cs.distance <= 4.5
                                          AND cs.is_same_entity = false
                                          AND h.ligand_id = $1
                                          AND cs.biomolecule_id = $2
                                          -- INTERACTION MUST BE WITH POLYMER ATOM (PROT/DNA/RNA)
                                          AND cs.structural_interaction_type_bm & 56 > 0
                                    UNION ALL
                                   SELECT DISTINCT h.ligand_id, cs.atom_bgn_id as atom_id
                                     FROM credo.hetatms h
                                     JOIN credo.contacts cs ON cs.atom_end_id = h.atom_id
                                    WHERE cs.distance <= 4.5
                                          AND cs.is_same_entity = false
                                          AND h.ligand_id = $1
                                          AND cs.biomolecule_id = $2
                                          -- INTERACTION MUST BE WITH POLYMER ATOM (PROT/DNA/RNA)
                                          AND cs.structural_interaction_type_bm & 3584 > 0
                                  ) sq
                         GROUP BY sq.ligand_id
                        )
                        SELECT fc.ligand_fragment_id, fc.num_int_atoms, f.num_hvy_atoms, W.num_int_atoms, l.num_hvy_atoms,
                               (fc.num_int_atoms / f.num_hvy_atoms::float) / (w.num_int_atoms / l.num_hvy_atoms::float) as fcd
                          FROM W
                          JOIN credo.ligand_fragments lf ON W.ligand_id = lf.ligand_id
                          JOIN credo.ligands l ON W.ligand_id = l.ligand_id
                          JOIN pdbchem.fragments f ON f.fragment_id = lf.fragment_id
                          JOIN credo.ligand_fcd fc ON lf.ligand_fragment_id = fc.ligand_fragment_id
                       ) sq
                 WHERE sq.ligand_fragment_id = f.ligand_fragment_id;
                ' USING lig_id, biomol_id;

                RAISE NOTICE 'updated fragment contact densities for all fragments of ligand %', lig_id;
        END LOOP;
END$$;

-- UPDATE THE FCD VALUE IN THE LIGAND_FRAGMENTS TABLE
UPDATE credo.ligand_fragments lf
   SET fcd = fcd.fcd
  FROM credo.ligand_fcd fcd
 WHERE lf.ligand_fragment_id = fcd.ligand_fragment_id;

-- SET LIGAND FRAGMENT IS_INTERACTING FLAG:
-- ONLY THOSE THAT HAVE FCD ARE INTERACTING
UPDATE credo.ligand_fragments lf
   SET is_interacting = true
  FROM credo.ligand_fcd fc
 WHERE fc.ligand_fragment_id = lf.ligand_fragment_id