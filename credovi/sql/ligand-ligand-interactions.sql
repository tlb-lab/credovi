  INSERT INTO credo.lig_lig_interactions(biomolecule_id, lig_bgn_id, lig_end_id)
  SELECT DISTINCT rp.biomolecule_id,
         lc1.ligand_id AS lig_bgn_id, lc2.ligand_id AS lig_end_id
    FROM credo.residue_interaction_pairs rp
    JOIN credo.ligand_components lc1 ON lc1.residue_id = rp.residue_bgn_id
    JOIN credo.ligand_components lc2 ON lc2.residue_id = rp.residue_end_id
   WHERE lc1.ligand_id != lc2.ligand_id
ORDER BY 1,2,3;

---- insert new ligand-ligand interactions
--DO $$
--    DECLARE
--        biomol_id INTEGER;
--    BEGIN
--        -- loop through all biomolecules that are not in the ligand-ligand interactions table yet
--        FOR biomol_id IN   SELECT biomolecule_id
--                             FROM credo.biomolecules
--                            WHERE biomolecule_id > (SELECT COALESCE(max(biomolecule_id),0)
--                                                      FROM credo.lig_lig_interactions)
--                         ORDER BY 1
--        LOOP
--            EXECUTE
--            '
--                 INSERT INTO credo.lig_lig_interactions(biomolecule_id, lig_bgn_id, lig_end_id, has_clash)
--                 SELECT DISTINCT cs.biomolecule_id,
--                        hbgn.ligand_id AS lig_bgn_id, hend.ligand_id AS lig_end_id,
--                        bool_or(cs.is_clash)
--                   FROM credo.contacts cs
--                   JOIN credo.hetatms hbgn ON hbgn.atom_id = cs.atom_bgn_id
--                   JOIN credo.hetatms hend ON hend.atom_id = cs.atom_end_id
--                  WHERE cs.biomolecule_id = $1
--                        -- no intramolecular contacts
--                        AND cs.is_same_entity = false
--                        AND hbgn.ligand_id != hend.ligand_id
--                        -- only feature contacts
--                        AND cs.distance <= 4.5
--               GROUP BY cs.biomolecule_id, hbgn.ligand_id, hend.ligand_id
--               ORDER BY 1,2,3;
--            ' USING biomol_id;
--
--            RAISE NOTICE 'inserted ligand-ligand interactions for biomolecule %', biomol_id;
--        END LOOP;
--END$$;

-- update ligand-ligand interactions that are not from the ASU
UPDATE credo.lig_lig_interactions li
   SET is_quaternary = true
  FROM credo.ligands l1, credo.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id AND l2.ligand_id = li.lig_end_id
       AND (l1.is_at_identity = false OR l2.is_at_identity = false);

-- update ligand-ligand interactions that have clashes
UPDATE credo.lig_lig_interactions li
   SET has_clash = true
  FROM credo.ligands l1, credo.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id AND l2.ligand_id = li.lig_end_id
       AND (l1.is_clashing = true OR l2.is_clashing = true);

-- update ligand-ligand interactions with identical ligands
UPDATE credo.lig_lig_interactions li
   SET is_homo_dimer = true
  FROM credo.ligands l1, credo.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND l1.ligand_name = l2.ligand_name;

-- update ligand-ligand interactions that involve an enzyme reaction product
UPDATE credo.lig_lig_interactions li
   SET has_product = true
  FROM credo.ligands l1, credo.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND (l1.is_product = true OR l2.is_product = true);

-- update ligand-ligand interactions with only drug-like ligands
UPDATE credo.lig_lig_interactions li
   SET has_drug_like_ligands = true
  FROM credo.ligands l1, credo.ligands l2,
       pdbchem.chem_comps c1, pdbchem.chem_comps c2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND l1.ligand_name = c1.het_id
       AND l2.ligand_name = c2.het_id
       AND c1.is_drug_like = true AND c2.is_drug_like = true;

-- update ligand-ligand interactions that involve an enzyme reaction substrate
UPDATE credo.lig_lig_interactions li
   SET has_substrate = true
  FROM credo.ligands l1, credo.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND (l1.is_substrate = true OR l2.is_substrate = true);

-- update ligand-ligand interactions with missing atoms
UPDATE credo.lig_lig_interactions li
   SET has_incomplete_res = true
  FROM credo.ligands l1, credo.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND (l1.is_incomplete = true OR l2.is_incomplete = true);

-- update ligand-ligand interactions that involve drug-target interactions
UPDATE credo.lig_lig_interactions li
   SET has_drug_target_int = true
  FROM credo.ligands l1, credo.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND (l1.is_drug_target_int = true OR l2.is_drug_target_int = true);

UPDATE credo.lig_lig_interactions li
   SET path = sq.path
  FROM (
        SELECT li.lig_lig_interaction_id,
                (
                 b.path || ('LLI:' || li.lig_lig_interaction_id) || lbgn.ligand_name || lend.ligand_name
                ) ::ptree AS path
          FROM credo.lig_lig_interactions li
          JOIN credo.biomolecules b USING(biomolecule_id)
          JOIN credo.ligands lbgn ON li.lig_bgn_id = lbgn.ligand_id
          JOIN credo.ligands lend ON li.lig_end_id = lend.ligand_id
       ) sq
 WHERE sq.lig_lig_interaction_id = li.lig_lig_interaction_id;
