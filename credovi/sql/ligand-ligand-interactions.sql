  INSERT INTO credo_dev.lig_lig_interactions(biomolecule_id, lig_bgn_id, lig_end_id)
  SELECT DISTINCT rp.biomolecule_id,
         lc1.ligand_id AS lig_bgn_id, lc2.ligand_id AS lig_end_id
    FROM credo_dev.residue_interaction_pairs rp
    JOIN credo_dev.ligand_components lc1 ON lc1.residue_id = rp.residue_bgn_id
    JOIN credo_dev.ligand_components lc2 ON lc2.residue_id = rp.residue_end_id
   WHERE lc1.ligand_id != lc2.ligand_id
ORDER BY 1,2,3;

-- update ligand-ligand interactions that are not from the ASU
UPDATE credo_dev.lig_lig_interactions li
   SET is_quaternary = true
  FROM credo_dev.ligands l1, credo_dev.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id AND l2.ligand_id = li.lig_end_id
       AND (l1.is_at_identity = false OR l2.is_at_identity = false);

-- update ligand-ligand interactions that have clashes
UPDATE credo_dev.lig_lig_interactions li
   SET has_clash = true
  FROM credo_dev.ligands l1, credo_dev.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id AND l2.ligand_id = li.lig_end_id
       AND (l1.is_clashing = true OR l2.is_clashing = true);

-- update ligand-ligand interactions with identical ligands
UPDATE credo_dev.lig_lig_interactions li
   SET is_homo_dimer = true
  FROM credo_dev.ligands l1, credo_dev.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND l1.ligand_name = l2.ligand_name;

-- update ligand-ligand interactions that involve an enzyme reaction product
UPDATE credo_dev.lig_lig_interactions li
   SET has_product = true
  FROM credo_dev.ligands l1, credo_dev.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND (l1.is_product = true OR l2.is_product = true);

-- update ligand-ligand interactions with only drug-like ligands
UPDATE credo_dev.lig_lig_interactions li
   SET has_drug_like_ligands = true
  FROM credo_dev.ligands l1, credo_dev.ligands l2,
       pdbchem.chem_comps c1, pdbchem.chem_comps c2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND l1.ligand_name = c1.het_id
       AND l2.ligand_name = c2.het_id
       AND c1.is_drug_like = true AND c2.is_drug_like = true;

-- update ligand-ligand interactions that involve an enzyme reaction substrate
UPDATE credo_dev.lig_lig_interactions li
   SET has_substrate = true
  FROM credo_dev.ligands l1, credo_dev.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND (l1.is_substrate = true OR l2.is_substrate = true);

-- update ligand-ligand interactions with missing atoms
UPDATE credo_dev.lig_lig_interactions li
   SET has_incomplete_res = true
  FROM credo_dev.ligands l1, credo_dev.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND (l1.is_incomplete = true OR l2.is_incomplete = true);

-- update ligand-ligand interactions that involve drug-target interactions
UPDATE credo_dev.lig_lig_interactions li
   SET has_drug_target_int = true
  FROM credo_dev.ligands l1, credo_dev.ligands l2
 WHERE l1.ligand_id = li.lig_bgn_id
       AND l2.ligand_id = li.lig_end_id
       AND (l1.is_drug_target_int = true OR l2.is_drug_target_int = true);

UPDATE credo_dev.lig_lig_interactions li
   SET path = sq.path
  FROM (
        SELECT li.lig_lig_interaction_id,
                (
                 b.path || ('LLI:' || li.lig_lig_interaction_id) || lbgn.ligand_name || lend.ligand_name
                ) ::ptree AS path
          FROM credo_dev.lig_lig_interactions li
          JOIN credo_dev.biomolecules b USING(biomolecule_id)
          JOIN credo_dev.ligands lbgn ON li.lig_bgn_id = lbgn.ligand_id
          JOIN credo_dev.ligands lend ON li.lig_end_id = lend.ligand_id
       ) sq
 WHERE sq.lig_lig_interaction_id = li.lig_lig_interaction_id;
