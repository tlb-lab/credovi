
    INSERT INTO credo_dev.lig_nuc_interactions(biomolecule_id, ligand_id, chain_nuc_id)
    SELECT DISTINCT rip.biomolecule_id, lc.ligand_id, n.chain_id
      FROM credo_dev.residue_interaction_pairs rip
      JOIN credo_dev.nucleotides n ON n.residue_id = rip.residue_bgn_id
      JOIN credo_dev.ligand_components lc ON lc.residue_id = rip.residue_end_id
     UNION
    SELECT DISTINCT rip.biomolecule_id, lc.ligand_id, n.chain_id
      FROM credo_dev.residue_interaction_pairs rip
      JOIN credo_dev.nucleotides n ON n.residue_id = rip.residue_end_id
      JOIN credo_dev.ligand_components lc ON lc.residue_id = rip.residue_bgn_id
  ORDER BY 1,2,3;


--
UPDATE credo_dev.lig_nuc_interactions li
   SET path = sq.path
  FROM (
        SELECT li.lig_nuc_interaction_id,
                (
                 b.path || ('LNI:' || li.lig_nuc_interaction_id) || l.ligand_name || c.pdb_chain_id
                ) ::ptree AS path
          FROM credo_dev.lig_nuc_interactions li
          JOIN credo_dev.biomolecules b USING(biomolecule_id)
          JOIN credo_dev.ligands l ON li.ligand_id = l.ligand_id
          JOIN credo_dev.chains c ON li.chain_nuc_id = c.chain_id
       ) sq
 WHERE sq.lig_nuc_interaction_id = li.lig_nuc_interaction_id;
