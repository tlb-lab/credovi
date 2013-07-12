-- CREATE A MAPPING BETWEEN LIGANDS AND FRAGMENTS
DO $$
    DECLARE
        lig_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL LIGANDS THAT DO NOT HAVE CLASHES AND ARE NOT DISORDERED
        FOR lig_id, biomol_id IN   SELECT ligand_id, biomolecule_id
                                     FROM credo_dev.ligands
                                    WHERE is_disordered = false
                                          -- no peptide ligands
                                          AND res_num IS NOT NULL
                                 ORDER BY 1
        LOOP
            EXECUTE
            '
            WITH ligand_fragments AS
                 (
                     INSERT INTO credo_dev.ligand_fragments(biomolecule_id, ligand_id, ligand_component_id, fragment_id, hit)
                     SELECT DISTINCT r.biomolecule_id, lc.ligand_id, lc.ligand_component_id, fragment_id, hit
                       FROM credo_dev.ligands l
                       JOIN credo_dev.ligand_components lc ON l.ligand_id = lc.ligand_id
                       JOIN credo_dev.residues r ON r.residue_id = lc.residue_id
                       JOIN pdbchem.chem_comp_fragments cf ON r.res_name = cf.het_id
                       JOIN pdbchem.chem_comp_fragment_atoms cfa ON cfa.chem_comp_fragment_id = cf.chem_comp_fragment_id
                      WHERE lc.ligand_id = $1
                   ORDER BY 1,2,3
                  RETURNING credo_dev.ligand_fragments.*
                 ),
                 ligand_fragment_atoms AS
                 (
                     INSERT INTO credo_dev.ligand_fragment_atoms(ligand_id, ligand_fragment_id, atom_id)
                     SELECT DISTINCT lf.ligand_id, lf.ligand_fragment_id, a.atom_id
                       FROM ligand_fragments lf
                            -- get the ligand component atoms
                       JOIN credo_dev.ligand_components lc ON lc.ligand_component_id = lf.ligand_component_id
                       JOIN credo_dev.residues r ON r.residue_id = lc.residue_id
                       JOIN credo_dev.atoms a on a.residue_id = lc.residue_id
                            -- link them to the fragment atoms
                       JOIN pdbchem.chem_comp_fragments cf
                            ON cf.het_id = r.res_name AND cf.fragment_id = lf.fragment_id
                       JOIN pdbchem.chem_comp_fragment_atoms cfa
                            ON cfa.chem_comp_fragment_id = cf.chem_comp_fragment_id
                            AND lf.hit = cfa.hit
                            AND cfa.pdb_name = a.atom_name
                      WHERE lc.ligand_id = lf.ligand_id AND a.biomolecule_id = $2
                   ORDER BY 1,2
                  RETURNING credo_dev.ligand_fragment_atoms.*
                 )
                 SELECT NULL;
            ' USING lig_id, biomol_id;
            RAISE NOTICE 'inserted ligand fragments and their atoms for ligand %', lig_id;
        END LOOP;
END$$;


--CREATE INDEX idx_ligand_fragments_biomolecule_id
--  ON credo_dev.ligand_fragments
--  USING btree
--  (biomolecule_id)
--  WITH (FILLFACTOR=100);

--CREATE INDEX idx_ligand_fragments_fragment_id
--  ON credo_dev.ligand_fragments
--  USING btree
--  (fragment_id, hit)
--  WITH (FILLFACTOR=100);

--CREATE INDEX idx_ligand_fragments_ligand_component_id
--  ON credo_dev.ligand_fragments
--  USING btree
--  (ligand_component_id)
--  WITH (FILLFACTOR=100);

--CREATE INDEX idx_ligand_fragments_ligand_id
--  ON credo_dev.ligand_fragments
--  USING btree
--  (ligand_id)
--  WITH (FILLFACTOR=100);

--CREATE INDEX idx_ligand_fragment_atoms_atom_id
--  ON credo_dev.ligand_fragment_atoms
--  USING btree
--  (atom_id)
--  WITH (FILLFACTOR=100);

--CREATE UNIQUE INDEX idx_ligand_fragment_atoms_ligand_fragment_id
--  ON credo_dev.ligand_fragment_atoms
--  USING btree
--  (ligand_fragment_id, atom_id)
--  WITH (FILLFACTOR=100);
--

--CREATE INDEX idx_ligand_fragment_atoms_ligand_id
--  ON credo_dev.ligand_fragment_atoms
--  USING btree
--  (ligand_id)
--  WITH (FILLFACTOR=100);


-- SET A FLAG OF ROOT LIGAND FRAGMENTS
UPDATE credo_dev.ligand_fragments lf
   SET is_root = true
  FROM credo_dev.ligand_components lc,
       pdbchem.fragment_hierarchies fh
 WHERE lc.ligand_component_id = lf.ligand_component_id
       AND fh.parent_id = lf.fragment_id
       AND fh.order_parent = 0
       AND fh.het_id = lc.het_id;

DO $$
    DECLARE
        lig_frag_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR lig_frag_id, biomol_id IN   SELECT lf.ligand_fragment_id, lf.biomolecule_id
                                          FROM credo_dev.ligand_fragments lf
                                          JOIN credo_dev.ligands USING(ligand_id)
                                      ORDER BY 1
        LOOP
            EXECUTE
                '
                   WITH fragment_contacts AS
                        (
                           SELECT fc.ligand_fragment_id,
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
                             FROM (
                                      SELECT lfa.ligand_fragment_id,
                                             cs.atom_bgn_id as atom_fragment_id, cs.atom_end_id as atom_id, cs.*
                                        FROM credo_dev.ligand_fragment_atoms lfa
                                             -- FRAGMENT ATOM IS BGN
                                        JOIN credo_dev.contacts cs ON cs.atom_bgn_id = lfa.atom_id
                                       WHERE lfa.ligand_fragment_id = $1
                                             AND cs.biomolecule_id = $2
                                             AND cs.is_same_entity = false
                                             AND cs.distance <= 4.5
                                             -- INTERACTION MUST BE WITH POLYMER ATOM (PROT/DNA/RNA)
                                             AND cs.structural_interaction_type_bm & 56 > 0
                                   UNION ALL
                                      SELECT lfa.ligand_fragment_id,
                                             cs.atom_end_id as atom_fragment_id, cs.atom_bgn_id as atom_id, cs.*
                                        FROM credo_dev.ligand_fragment_atoms lfa
                                             -- FRAGMENT ATOM IS END
                                        JOIN credo_dev.contacts cs ON cs.atom_end_id = lfa.atom_id
                                       WHERE lfa.ligand_fragment_id = $1
                                             AND cs.biomolecule_id = $2
                                             AND cs.is_same_entity = false
                                             AND cs.distance <= 4.5
                                             -- INTERACTION MUST BE WITH POLYMER ATOM (PROT/DNA/RNA)
                                             AND cs.structural_interaction_type_bm & 3584 > 0
                                  ) fc
                         GROUP BY fc.ligand_fragment_id
                         ORDER BY fc.ligand_fragment_id
                        )
                 UPDATE credo_dev.ligand_fragments lf
                    SET num_int_atoms = fc.num_int_atoms,
                        num_contacts = fc.num_covalent + fc.num_vdw_clash + fc.num_vdw + fc.num_proximal,
                        num_covalent = fc.num_covalent,
                        num_vdw_clash = fc.num_vdw_clash,
                        num_vdw = fc.num_vdw,
                        num_proximal = fc.num_proximal,
                        num_hbond = fc.num_hbond,
                        num_weak_hbond = fc.num_weak_hbond,
                        num_xbond = fc.num_xbond,
                        num_ionic = fc.num_ionic,
                        num_metal_complex = fc.num_metal_complex,
                        num_aromatic = fc.num_aromatic,
                        num_hydrophobic = fc.num_hydrophobic,
                        num_carbonyl = fc.num_carbonyl
                   FROM fragment_contacts fc
                  WHERE fc.ligand_fragment_id = lf.ligand_fragment_id;
                  ' USING lig_frag_id, biomol_id;

                RAISE NOTICE 'updated fragment contact densities for ligand fragment %', lig_frag_id;
        END LOOP;
END$$;

UPDATE credo_dev.ligand_fragments u
   SET fcd = sq.fcd, fad = sq.fad
  FROM (
        -- get the number of ligand heavy atoms and interacting atoms from the root fragments
        WITH ligands AS
        (
         SELECT ligand_id,
                num_hvy_atoms, num_int_atoms, num_contacts,
                num_int_atoms / num_hvy_atoms::numeric AS int_ratio,
                num_contacts / num_hvy_atoms::numeric AS cs_ratio
           FROM credo_dev.ligand_fragments lf
           JOIN credo_dev.ligands l USING(ligand_id)
          WHERE lf.is_root = true and l.num_hvy_atoms > 0
        )
        SELECT lf.ligand_fragment_id,
               (lf.num_int_atoms / f.num_hvy_atoms::numeric) / l.int_ratio as fad,
               (lf.num_contacts / f.num_hvy_atoms::numeric) / l.cs_ratio as fcd
          FROM credo_dev.ligand_fragments lf
          JOIN ligands l ON l.ligand_id = lf.ligand_id
          JOIN pdbchem.fragments f ON f.fragment_id = lf.fragment_id
          WHERE f.num_hvy_atoms > 0 and l.int_ratio > 0 and l.cs_ratio > 0
       ) sq
 WHERE u.ligand_fragment_id = sq.ligand_fragment_id;

DO $$
    DECLARE
        lig_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        FOR lig_id, biomol_id IN   SELECT DISTINCT ligand_id, biomolecule_id
                                     FROM credo_dev.ligand_fragments
                                 ORDER BY 1
        LOOP
            EXECUTE
            '
            UPDATE credo_dev.ligand_fragments lf
               SET npr1 = (sq.nprs).npr1, npr2 = (sq.nprs).npr2
              FROM (
                      WITH sq AS
                           (
                               SELECT lfa.ligand_fragment_id,
                                      array_cat(a.coords) AS coords,
                                      array_agg(a.element) AS elements,
                                      count(a.atom_id) AS num_atoms
                                 FROM credo_dev.ligand_fragment_atoms lfa
                                 JOIN credo_dev.atoms a ON a.atom_id = lfa.atom_id
                                WHERE ligand_id = $1 AND a.biomolecule_id = $2
                             GROUP BY lfa.ligand_fragment_id
                           )
                    SELECT ligand_fragment_id, openeye.nprs(sq.coords, sq.elements) as nprs
                      FROM sq
                     WHERE sq.num_atoms >= 5
                   ) sq
             WHERE lf.ligand_fragment_id = sq.ligand_fragment_id;
            ' USING lig_id, biomol_id;
            RAISE NOTICE 'updated normalized ratios of PMI for ligand %', lig_id;
        END LOOP;
END$$;

DO $$
    DECLARE
        lig_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        FOR lig_id, biomol_id IN   SELECT DISTINCT ligand_id, biomolecule_id
                                     FROM credo_dev.ligand_fragments
                                 ORDER BY 1
        LOOP
            EXECUTE
            '
                 WITH atoms AS
                      (
                       SELECT ligand_fragment_id, atom_id,
                              asa_apo, asa_bound, asa_delta,
                              CASE WHEN a.element = ''N''
                                        OR a.element = ''O''
                                        OR a.element = ''P''
                                        OR a.element = ''S''
                                   THEN true ELSE false
                              END AS is_polar
                         FROM credo_dev.binding_site_atom_surface_areas bsa
                         JOIN credo_dev.ligand_fragment_atoms lfa USING(ligand_id, atom_id)
                         JOIN credo_dev.atoms a USING(atom_id)
                        WHERE ligand_id = $1 AND a.biomolecule_id = $2
                      ),
                      asa AS
                      (
                         SELECT ligand_fragment_id,
                                sum(asa_apo) as asa_apo,
                                sum(asa_bound) as asa_bound,
                                sum(asa_delta) as asa_delta,
                                -- catch division by zero error if all fragment
                                -- atoms are already buried
                                CASE WHEN sum(asa_apo) > 0
                                     THEN sum(asa_delta) / sum(asa_apo)
                                     ELSE NULL
                                END AS asa_buriedness
                           FROM atoms
                       GROUP BY ligand_fragment_id
                      ),
                      aasa AS
                      (
                         SELECT ligand_fragment_id,
                                sum(asa_apo) as aasa_apo,
                                sum(asa_bound) as aasa_bound,
                                sum(asa_delta) as aasa_delta,
                                CASE WHEN sum(asa_apo) > 0
                                     THEN sum(asa_delta) / sum(asa_apo)
                                     ELSE NULL
                                END AS aasa_buriedness
                           FROM atoms
                          WHERE is_polar = false
                       GROUP BY ligand_fragment_id
                      ),
                      pasa AS
                      (
                         SELECT ligand_fragment_id,
                                sum(asa_apo) as pasa_apo,
                                sum(asa_bound) as pasa_bound,
                                sum(asa_delta) as pasa_delta,
                                CASE WHEN sum(asa_apo) > 0
                                     THEN sum(asa_delta) / sum(asa_apo)
                                     ELSE NULL
                                END AS pasa_buriedness
                           FROM atoms
                          WHERE is_polar = true
                       GROUP BY ligand_fragment_id
                      )
               UPDATE credo_dev.ligand_fragments lf
                  SET asa_buriedness = sq.asa_buriedness,
                      aasa_buriedness = sq.aasa_buriedness,
                      pasa_buriedness = sq.pasa_buriedness
                 FROM (
                          SELECT asa.ligand_fragment_id,
                                 asa_buriedness, aasa_buriedness, pasa_buriedness
                            FROM asa
                       LEFT JOIN aasa ON asa.ligand_fragment_id = aasa.ligand_fragment_id
                       LEFT JOIN pasa ON asa.ligand_fragment_id = pasa.ligand_fragment_id
                      ) sq
                WHERE lf.ligand_fragment_id = sq.ligand_fragment_id
            ' USING lig_id, biomol_id;
            RAISE NOTICE 'Updated surface buriedness for all fragments of ligand %', lig_id;
        END LOOP;
END$$;

UPDATE credo_dev.ligand_fragments lf
   SET rel_asa_buriedness = sq.rel_asa_buriedness
  FROM (
          with roots AS
               (
                SELECT ligand_id, asa_buriedness
                  FROM credo_dev.ligand_fragments lf
                 WHERE is_root = true AND asa_buriedness > 0
               )
        SELECT ligand_fragment_id, lf.asa_buriedness / r.asa_buriedness as rel_asa_buriedness
          FROM credo_dev.ligand_fragments lf
          JOIN roots r ON r.ligand_id = lf.ligand_id
       ) sq
 WHERE lf.ligand_fragment_id = sq.ligand_fragment_id;
