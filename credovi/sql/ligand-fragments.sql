DROP    TABLE IF EXISTS credo.ligand_fragments;
CREATE  TABLE credo.ligand_fragments (
        ligand_fragment_id SERIAL NOT NULL,
        ligand_id INTEGER NOT NULL,
        ligand_component_id INTEGER NOT NULL,
        fragment_id INTEGER NOT NULL,
        hit INTEGER NOT NULL,
        usr_space CUBE DEFAULT NULL,
        usr_moments_env FLOAT[] DEFAULT NULL,
        is_root boolean NOT NULL DEFAULT FALSE,
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

CREATE  INDEX idx_ligand_fragments_usr_space
        ON credo.ligand_fragments
        USING gist (usr_space);

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

GRANT   SELECT ON TABLE credo.ligand_fcd TO credouser;

-- CREATE A MAPPING BETWEEN LIGANDS AND FRAGMENTS
INSERT      INTO credo.ligand_fragments(ligand_id, ligand_component_id, fragment_id, hit)
SELECT      DISTINCT lc.ligand_id, lc.ligand_component_id, fragment_id, hit
FROM        credo.ligand_components lc
JOIN        credo.residues r ON r.residue_id = lc.residue_id
JOIN        pdbchem.chem_comp_fragments cf ON r.res_name = cf.het_id
JOIN        pdbchem.chem_comp_fragment_atoms cfa ON cfa.chem_comp_fragment_id = cf.chem_comp_fragment_id
ORDER BY    1,2,3;


-- CREATE A MAPPING BETWEEN LIGAND FRAGMENTS AND THEIR ATOMS
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR biomol_id IN SELECT biomolecule_id FROM credo.biomolecules ORDER BY 1 
        LOOP
            -- INSERT LIGAND FRAGMENT ATOMS
            EXECUTE 
            '
            INSERT      INTO credo.ligand_fragment_atoms(ligand_fragment_id, atom_id)
            SELECT      DISTINCT lf.ligand_fragment_id, a.atom_id
            FROM        credo.ligand_fragments lf
                        -- GET THE LIGAND COMPONENT ATOMS
            JOIN        credo.ligand_components lc ON lc.ligand_component_id = lf.ligand_component_id
            JOIN        credo.residues r ON r.residue_id = lc.residue_id
            JOIN        credo.atoms a on a.residue_id = lc.residue_id
                        -- LINK THEM TO THE FRAGMENT ATOMS
            JOIN        pdbchem.chem_comp_fragments cf 
                        ON cf.het_id = r.res_name AND cf.fragment_id = lf.fragment_id
            JOIN        pdbchem.chem_comp_fragment_atoms cfa
                        ON cfa.chem_comp_fragment_id = cf.chem_comp_fragment_id
                        AND lf.hit = cfa.hit
                        AND cfa.pdb_name = a.atom_name
             WHERE      a.biomolecule_id = $1
            ORDER BY    1,2;
            ' USING biomol_id;

            RAISE NOTICE 'inserted ligand fragments for biomolecule %', biomol_id;
        END LOOP;
END$$;

-- FRAGMENT CONTACT DENSITY FOR EACH LIGAND FRAGMENT

DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR biomol_id IN SELECT biomolecule_id FROM credo.biomolecules ORDER BY 1 
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
                            SELECT      DISTINCT lfa.ligand_fragment_id,
                                        cs.atom_bgn_id as atom_fragment_id, cs.atom_end_id as atom_id, cs.*
                            FROM        credo.ligand_fragment_atoms lfa
                                        -- FRAGMENT ATOM IS BGN
                            JOIN        credo.contacts cs ON cs.atom_bgn_id = lfa.atom_id
                            WHERE       cs.biomolecule_id = $1
                                        AND cs.is_secondary = false 
                                        AND cs.is_same_entity = false 
                                        AND cs.distance <= 4.5
                            UNION ALL
                            SELECT      DISTINCT lfa.ligand_fragment_id,
                                        cs.atom_end_id as atom_fragment_id, cs.atom_bgn_id as atom_id, cs.*
                            FROM        credo.ligand_fragment_atoms lfa
                                        -- FRAGMENT ATOM IS END
                            JOIN        credo.contacts cs ON cs.atom_end_id = lfa.atom_id
                            WHERE       cs.biomolecule_id = $1
                                        AND cs.is_secondary = false 
                                        AND cs.is_same_entity = false 
                                        AND cs.distance <= 4.5
                            )
                SELECT      fc.ligand_fragment_id,
                            COUNT(DISTINCT fc.atom_id) AS num_int_atoms,
                            SUM(ARRAY_AGG(fc.is_covalent)) as num_covalent,
                            SUM(ARRAY_AGG(fc.is_vdw_clash)) as num_vdw_clash,
                            SUM(ARRAY_AGG(fc.is_vdw)) AS num_vdw,
                            SUM(ARRAY_AGG(fc.is_proximal)) AS num_proximal,
                            SUM(ARRAY_AGG(fc.is_hbond)) AS num_hbond,
                            SUM(ARRAY_AGG(fc.is_weak_hbond)) AS num_weak_hbond,
                            SUM(ARRAY_AGG(fc.is_xbond)) as num_xbond,
                            SUM(ARRAY_AGG(fc.is_ionic)) as num_ionic,
                            SUM(ARRAY_AGG(fc.is_metal_complex)) AS num_metal_complex,
                            SUM(ARRAY_AGG(fc.is_aromatic)) AS num_aromatic,
                            SUM(ARRAY_AGG(fc.is_hydrophobic)) AS num_hydrophobic,
                            SUM(ARRAY_AGG(fc.is_carbonyl)) AS num_carbonyl
                FROM        fragment_contacts fc
                GROUP BY    fc.ligand_fragment_id
                ORDER BY    fc.ligand_fragment_id
                ' USING biomol_id;
                
                RAISE NOTICE 'calculated fragment contact densities for biomolecule %', biomol_id;
        END LOOP;
END$$;

-- CALCULATE FCD
UPDATE credo.ligand_fcd f
   SET fcd = sq.fcd
  FROM (
        WITH W AS
        (
           SELECT sq.ligand_id, sq.num_hvy_atoms, COUNT(DISTINCT sq.atom_id) AS num_int_atoms
             FROM (
                   SELECT DISTINCT l.ligand_id, l.num_hvy_atoms, cs.atom_end_id as atom_id
                     FROM credo.ligand_fragments lf
                     JOIN credo.ligands l USING(ligand_id)
                     JOIN credo.hetatms h USING(ligand_id)
                     JOIN credo.contacts cs ON cs.atom_bgn_id = h.atom_id
                    WHERE l.is_incomplete = false
                          AND cs.is_secondary = false 
                          AND cs.distance <= 4.5

                    UNION
                   SELECT DISTINCT l.ligand_id, l.num_hvy_atoms, cs.atom_bgn_id
                     FROM credo.ligand_fragments lf
                     JOIN credo.ligands l USING(ligand_id)
                     JOIN credo.hetatms h USING(ligand_id)
                     JOIN credo.contacts cs ON cs.atom_end_id = h.atom_id
                    WHERE l.is_incomplete = false
                          AND cs.is_secondary = false AND cs.distance <= 4.5
                  ) sq
         GROUP BY sq.ligand_id, sq.num_hvy_atoms
        )
        SELECT      fc.ligand_fragment_id, fc.num_int_atoms, f.num_atoms, W.num_int_atoms, W.num_hvy_atoms,
                    (fc.num_int_atoms / f.num_atoms::float) / (w.num_int_atoms / W.num_hvy_atoms::float) as fcd
        FROM        W
        JOIN        credo.ligand_fragments lf ON W.ligand_id = lf.ligand_id
        JOIN        pdbchem.fragments f ON f.fragment_id = lf.fragment_id
        JOIN        credo.ligand_fcd fc ON lf.ligand_fragment_id = fc.ligand_fragment_id
        ) sq
 WHERE sq.ligand_fragment_id = f.ligand_fragment_id;

/*
-- UPDATE THE USR ENVIRONMENT MOMENTS FOR EACH LIGAND FRAGMENT
DROP        TABLE IF EXISTS tmp_ligand_fragment_env;
CREATE      TEMP TABLE tmp_ligand_fragment_env AS
SELECT      DISTINCT lfa.ligand_fragment_id, a2.coords as atoms,
            CASE WHEN a2.is_hydrophobe THEN a2.coords ELSE NULL END as hydrophobes,
            CASE WHEN a2.is_aromatic THEN a2.coords ELSE NULL END as aromatics,
            CASE WHEN a2.is_acceptor THEN a2.coords ELSE NULL END as acceptors,
            CASE WHEN a2.is_donor THEN a2.coords ELSE NULL END as donors
FROM        credo.ligand_fragment_atoms lfa
JOIN        credo.ligand_fragments USING(ligand_fragment_id)
JOIN        pdbchem.fragments f USING(fragment_id)
JOIN        credo.contacts cs ON cs.atom_end_id = lfa.atom_id
JOIN        credo.atoms a1 ON a1.atom_id = cs.atom_bgn_id
            -- GET ALL RESIDUE ATOMS
JOIN        credo.atoms a2 USING(residue_id)
WHERE       -- ONLY KEEP ENVIRONMENT MOMENTS FOR FRAGMENTS WITH AT LEAST 5 ATOMS
            f.num_atoms >= 5
UNION
SELECT      DISTINCT lfa.ligand_fragment_id, a2.coords as atoms,
            CASE WHEN a2.is_hydrophobe THEN a2.coords ELSE NULL END as hydrophobes,
            CASE WHEN a2.is_aromatic THEN a2.coords ELSE NULL END as aromatics,
            CASE WHEN a2.is_acceptor THEN a2.coords ELSE NULL END as acceptors,
            CASE WHEN a2.is_donor THEN a2.coords ELSE NULL END as donors
FROM        credo.ligand_fragment_atoms lfa
JOIN        credo.ligand_fragments USING(ligand_fragment_id)
JOIN        pdbchem.fragments f USING(fragment_id)
JOIN        credo.contacts cs ON cs.atom_bgn_id = lfa.atom_id
JOIN        credo.atoms a1 ON a1.atom_id = cs.atom_end_id
            -- GET ALL RESIDUE ATOMS
JOIN        credo.atoms a2 USING(residue_id)
WHERE       -- ONLY KEEP ENVIRONMENT MOMENTS FOR FRAGMENTS WITH AT LEAST 5 ATOMS
            f.num_atoms >= 5;

-- CREATE INDEX FOR THE TEMP TABLE
CREATE  INDEX idx_tmp_ligand_fragment_env ON tmp_ligand_fragment_env USING btree (ligand_fragment_id);

UPDATE      credo.ligand_fragments lf
SET         usr_moments_env = sq.usr_moments_env
FROM        (
            SELECT      ligand_fragment_id,
                        coords_to_usr(ARRAY_AGG(atoms) || ARRAY_AGG(hydrophobes) || ARRAY_AGG(aromatics) ||
                                      ARRAY_AGG(acceptors) || ARRAY_AGG(donors)) as usr_moments_env
            FROM        tmp_ligand_fragment_env
            GROUP BY    ligand_fragment_id
                        -- ONLY KEEP ENVIRONMENTS WITH AT LEAST 10 ATOMS
            HAVING      ARRAY_LENGTH(ARRAY_AGG(atoms),1) >= 10
            ) sq
WHERE       lf.ligand_fragment_id = sq.ligand_fragment_id;

-- UPDATE USR SPACE CUBE FOR USR ENVIRONMENT MOMENTS
UPDATE      credo.ligand_fragments lf
SET         usr_space = cube(lf.usr_moments_env[1:12]);
*/

