-- STRUCTURE_ID IN MMCIF IS PDB CODE!
-- MMCIF.pdb_strand_to_entity_id IS REQUIRED!
-- MAKE SURE TO JOIN ON chain.pdb_chain_asu!

-- UPDATE THE PDB STRUCTURE TITLE
UPDATE      credo.structures s
SET         title = t.title
FROM        mmcif.struct t
WHERE       s.pdb = t.structure_id;

-- UPDATE AUTHORS OF PDB ENTRY
UPDATE      credo.structures s
SET         authors = rs.authors
FROM        (
            select      structure_id, string_agg(name,', ') as authors
            from        mmcif.citation_author
            where       citation_id = 'primary'
            group by    structure_id
            ) rs
WHERE       s.pdb = rs.structure_id;

-- UPDATE DEPOSITION AND LATEST MOD DATE
UPDATE credo.structures s
   SET deposition_date = sq.deposition_date, modified_date = sq.modified_date
  FROM (
          SELECT structure_id as pdb, min(date_original) as deposition_date, max(date) as modified_date
            FROM mmcif.database_pdb_rev
        GROUP BY structure_id
       ) sq
 WHERE s.pdb = sq.pdb;

-- UPDATE THE EXPERIMENTAL METHOD FOR EACH STRUCTURE
UPDATE      credo.structures s
SET         method = e.method
FROM        mmcif.exptl e
WHERE       s.pdb = e.structure_id;

-- UPDATE IMPORTANT STRUCTURE FACTORS
UPDATE      credo.structures s
SET         resolution = rs.resolution,
            r_factor = rs.r_factor,
            r_free = rs.r_free
FROM        (
            SELECT  s.structure_id,
                    r.ls_d_res_high as resolution,
                    NULLIF(r.ls_R_factor_R_work, 0) as r_factor,
                    NULLIF(r.ls_R_factor_R_free, 0) as r_free
            FROM    credo.structures s
            JOIN    mmcif.refine r ON r.structure_id = s.pdb
            ) rs
WHERE       s.structure_id = rs.structure_id;

-- UPDATE PH OF CRYSTALLISATION
UPDATE      credo.structures s
SET         ph = c.ph
FROM        mmcif.exptl_crystal_grow c
WHERE       s.pdb = c.structure_id;

-- UPDATE STRUCTURE DIFFRACTION-COMPONENT-PRECISION INDEX (DPI)

UPDATE      credo.structures s
SET         dpi = rs.dpi,
            dpi_theoretical_min = rs.dpi_theoretical_min
FROM        (
            SELECT 	s.structure_id,
                    2.2 * SQRT(x1.natoms) * CBRT(x2.va) * POW(x3.nobs, -5/6.0) * x3.r_free AS dpi,
                    0.22 * POW(2.4, -1/2.0) * POW(1.0, -5/6.0) * x3.r_free * POW(x3.resolution, 5/2.0) AS dpi_theoretical_min
            FROM 	credo.structures s,
                    -- NUMBER OF ATOMS
                    (
                    SELECT 		structure_id, number_atoms_total AS natoms
                    FROM 		mmcif.refine_hist
                    ) x1,
                    -- UNIT CELL VOLUME
                    (
                    SELECT 		structure_id, (length_a * length_b * length_c) AS va
                    FROM 		mmcif.cell1
                    ) x2,
                    -- R FREE AND NUMBER OF CRYSTALLOGRAPHIC OBSERVATIONS
                    (
                    SELECT 		structure_id, ls_number_reflns_obs AS nobs,
                                ls_R_factor_R_free AS r_free,
                                ls_d_res_high AS resolution
                    FROM 		mmcif.refine
                    ) x3
            WHERE 	s.pdb = x1.structure_id
                    AND x1.structure_id = x2.structure_id
                    AND x2.structure_id = x3.structure_id
                    AND x3.r_free > 0 AND x3.nobs > 0
                    AND x3.r_free BETWEEN 0 AND 1
            ) rs
WHERE       s.structure_id = rs.structure_id;

-- SET THE ENTITY TYPE FOR EACH CHAIN
UPDATE      credo.chains c
SET         chain_type = rs.type
FROM        (
            SELECT  DISTINCT c.chain_id, ep.type
            FROM    mmcif.entity_poly ep
            JOIN    mmcif.pdb_strand_to_entity_id pse
                    ON pse.structure_id = ep.structure_id AND pse.entity_id = ep.entity_id
            JOIN    credo.structures s ON s.pdb = pse.structure_id
            JOIN    credo.biomolecules b ON b.structure_id = s.structure_id
            JOIN    credo.chains c
                    ON c.biomolecule_id = b.biomolecule_id
                    AND c.pdb_chain_asu_id = pse.pdb_strand_id
            ) rs
WHERE       c.chain_id = rs.chain_id;

-- UPDATE TITLE OF PDB CHAINS
UPDATE      credo.chains c
SET         title = TRIM(rs.pdbx_description)
FROM        (
            SELECT  DISTINCT c.chain_id, e.pdbx_description
            FROM    mmcif.entity e
            JOIN    mmcif.pdb_strand_to_entity_id pse
                    ON pse.structure_id = e.structure_id AND pse.entity_id = e.id
            JOIN    credo.structures s ON s.pdb = pse.structure_id
            JOIN    credo.biomolecules b ON b.structure_id = s.structure_id
            JOIN    credo.chains c
                    ON c.biomolecule_id = b.biomolecule_id
                    AND c.pdb_chain_asu_id = pse.pdb_strand_id
            ) rs
WHERE       c.chain_id = rs.chain_id;

-- UPDATE PDB CHAIN POLYMER SEQUENCE
UPDATE      credo.chains c
SET         chain_seq = rs.pdbx_seq_one_letter_code_can,
            chain_seq_md5 = UPPER(md5(rs.pdbx_seq_one_letter_code_can)),
            chain_length = LENGTH(rs.pdbx_seq_one_letter_code_can)
FROM        (
            SELECT  DISTINCT c.chain_id, p.pdbx_seq_one_letter_code_can
            FROM    mmcif.entity_poly p
            JOIN    mmcif.pdb_strand_to_entity_id pse
                    ON pse.structure_id = p.structure_id AND pse.entity_id = p.entity_id
            JOIN    credo.structures s ON s.pdb = pse.structure_id
            JOIN    credo.biomolecules b ON b.structure_id = s.structure_id
            JOIN    credo.chains c
                    ON c.biomolecule_id = b.biomolecule_id
                    AND c.pdb_chain_asu_id = pse.pdb_strand_id
            ) rs
WHERE       c.chain_id = rs.chain_id;

-- INSERT POLYPEPTIDES
   INSERT INTO credo.polypeptides(chain_id)
   SELECT c.chain_id
     FROM credo.chains c
          -- IGNORE ALREADY EXISTING ONES
    WHERE c.chain_id > (SELECT COALESCE(MAX(chain_id),0) FROM credo.polypeptides)
          AND (c.chain_type = 'polypeptide(L)' OR c.chain_type = 'polypeptide(D)')
 ORDER BY 1;

-- ENZYME FLAG
UPDATE  credo.polypeptides p
SET     is_enzyme = true
FROM    (
        SELECT  DISTINCT c.chain_id
        FROM    credo.chains c
        JOIN    credo.biomolecules b USING(biomolecule_id)
        JOIN    credo.structures s USING(structure_id)
        JOIN    pdb.map_regions m ON s.pdb = m.pdb AND c.pdb_chain_asu_id = m.pdb_chain_id
        WHERE   m.db_source = 'EC'
        ) sq
WHERE   p.chain_id = sq.chain_id;

-- HUMAN FLAG
UPDATE  credo.polypeptides p
SET     is_human = true
FROM    (
        SELECT  DISTINCT c.chain_id
        FROM    credo.chains c
        JOIN    credo.biomolecules b USING(biomolecule_id)
        JOIN    credo.structures s USING(structure_id)
        JOIN    pdb.map_regions m ON s.pdb = m.pdb AND c.pdb_chain_asu_id = m.pdb_chain_id
        WHERE   m.db_source = 'NCBI' AND m.db_accession_id = '9606'
        ) sq
WHERE   p.chain_id = sq.chain_id;

-- DRUG TARGET FLAG
UPDATE  credo.polypeptides p
SET     is_drug_target = true
FROM    (
        SELECT  DISTINCT c.chain_id
        FROM    credo.chains c
        JOIN    credo.biomolecules b USING(biomolecule_id)
        JOIN    credo.structures s USING(structure_id)
        JOIN    pdb.map_regions m ON s.pdb = m.pdb AND c.pdb_chain_asu_id = m.pdb_chain_id
        JOIN    chembl.target_dictionary td ON td.protein_accession = m.db_accession_id
        WHERE   m.db_source = 'UniProt'
        ) sq
WHERE   p.chain_id = sq.chain_id;

-- SET A FLAG FOR POLYPEPTIDE THAT ARE IN UNIPROT
  WITH chains AS
       (
        SELECT DISTINCT p.chain_id
          FROM credo.polypeptides p
          JOIN credo.xrefs xr ON xr.entity_type = 'Chain' AND p.chain_id = xr.entity_id
         WHERE xr.source = 'UniProt'
       )
UPDATE credo.polypeptides p
   SET is_in_uniprot = true
  FROM chains c
 WHERE c.chain_id = p.chain_id;

-- set a flag for polypeptide chains that are kinases
  WITH uniprots AS
       (
        SELECT acc_human as uniprot FROM credo.kinases WHERE acc_human IS NOT NULL
        UNION
        SELECT acc_mouse FROM credo.kinases WHERE acc_mouse IS NOT NULL
       ),
       chains AS
       (
        SELECT DISTINCT p.chain_id
          FROM credo.polypeptides p
          JOIN credo.xrefs xr ON xr.entity_type = 'Chain' and p.chain_id = xr.entity_id
          JOIN uniprots u ON xr.source = 'UniProt' AND u.uniprot = xr.xref
       )
UPDATE credo.polypeptides p
   SET is_kinase = true
  FROM chains c
 WHERE c.chain_id = p.chain_id;

-- INSERT OLIGONUCLEOTIDES
INSERT      INTO credo.oligonucleotides(chain_id, nucleic_acid_type)
SELECT      c.chain_id,
            CASE
                WHEN c.chain_type = 'polyribonucleotide' THEN 'RNA'
                WHEN c.chain_type = 'polydeoxyribonucleotide' THEN 'DNA'
                WHEN c.chain_type = 'polydeoxyribonucleotide/polyribonucleotide hybrid' THEN 'HYB'
            END as nucleic_acid_type
FROM        credo.chains c
LEFT JOIN   credo.oligonucleotides n ON n.chain_id = c.chain_id
WHERE       n.chain_id IS NULL
            AND (c.chain_type = 'polyribonucleotide'
            OR c.chain_type = 'polydeoxyribonucleotide'
            OR c.chain_type = 'polydeoxyribonucleotide/polyribonucleotide hybrid')
ORDER BY    1;

-- INSERT POLYSACCHARIDES
INSERT      INTO credo.polysaccharides(chain_id)
SELECT      c.chain_id
FROM        credo.chains c
LEFT JOIN   credo.polysaccharides p ON p.chain_id = c.chain_id
WHERE       p.chain_id IS NULL AND c.chain_type = 'polysaccharide(D)'
ORDER BY    1;

-- UPDATE RESIDUE ONE-LETTER CODE
UPDATE  credo.peptides p
SET     one_letter_code = openeye.one_letter_code(res_name)
WHERE   one_letter_code IS NULL;

-- UPDATE NON-STANDARD RESIDUES
UPDATE  credo.peptides p
SET     is_non_std = true
WHERE   p.res_name not in ('CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 'THR', 'PHE', 'ALA',
                           'GLY', 'HIS', 'GLU', 'LEU', 'ARG', 'TRP', 'VAL', 'ASN', 'TYR', 'MET');

-- UPDATE RESIDUE MAP ID
UPDATE  credo.peptides p
SET     res_map_id = m.res_map_id
FROM    credo.chains c,
        credo.biomolecules b,
        credo.structures s,
        pdb.res_map m
WHERE   p.chain_id = c.chain_id
        AND c.biomolecule_id = b.biomolecule_id
        AND b.structure_id = s.structure_id
        AND m.pdb = s.pdb
        AND m.pdb_chain_id = c.pdb_chain_asu_id
        AND m.pdb_res_num = p.res_num
        AND m.pdb_ins_code = p.ins_code;

-- UPDATE MODIFIED RESIDUES
UPDATE  credo.peptides p
SET     is_modified = true
FROM    pdb.res_map m
WHERE   m.res_map_id = p.res_map_id
        AND m.is_modified = true;

-- MUTATED PEPTIDES
UPDATE 	credo.peptides p
SET     is_mutated = true
FROM    pdb.res_map m
WHERE   m.res_map_id = p.res_map_id
        AND p.is_modified = false AND p.is_non_std = false
        AND p.one_letter_code != m.uniprot_res_name;

-- UPDATE CATH FOR PEPTIDES
UPDATE  credo.peptides p
SET     cath = mr.db_accession_id
FROM    pdb.res_map m, pdb.map_regions mr
WHERE   m.res_map_id = p.res_map_id
        AND m.pdb = mr.pdb
        AND m.pdb_chain_id = mr.pdb_chain_id
        AND m.sifts_res_num BETWEEN mr.sifts_res_num_start AND mr.sifts_res_num_end
        AND mr.db_source = 'CATH';

-- UPDATE SCOP PX FOR PEPTIDES
UPDATE  credo.peptides p
SET     px = mr.db_accession_id::integer
FROM    pdb.res_map m, pdb.map_regions mr
WHERE   m.res_map_id = p.res_map_id
        AND m.pdb = mr.pdb
        AND m.pdb_chain_id = mr.pdb_chain_id
        AND m.sifts_res_num BETWEEN mr.sifts_res_num_start AND mr.sifts_res_num_end
        AND mr.db_source = 'SCOP';

-- PROTEIN-PROTEIN INTERFACES
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES THAT ARE NOT IN THE INTERFACES TABLE YET
        FOR biomol_id IN   SELECT biomolecule_id
                             FROM credo.biomolecules
                            WHERE biomolecule_id > (SELECT COALESCE(max(biomolecule_id),0) FROM credo.interfaces)
                         ORDER BY 1
        LOOP
            EXECUTE
            '
              INSERT INTO credo.interfaces(biomolecule_id, chain_bgn_id, chain_end_id,
                                           num_res_bgn, num_res_end)
              SELECT cs.biomolecule_id,
                     rbgn.chain_id as chain_bgn_id, rend.chain_id as chain_end_id,
                     COUNT(DISTINCT rbgn.residue_id) AS num_res_bgn,
                     COUNT(DISTINCT rend.residue_id) AS num_res_end
                FROM credo.contacts cs
                JOIN credo.atoms abgn
                     ON abgn.atom_id = cs.atom_bgn_id
                     AND abgn.biomolecule_id = cs.biomolecule_id
                JOIN credo.atoms aend
                     ON aend.atom_id = cs.atom_end_id
                     AND aend.biomolecule_id = cs.biomolecule_id
                JOIN credo.peptides rbgn ON rbgn.residue_id = abgn.residue_id
                JOIN credo.peptides rend ON rend.residue_id = aend.residue_id
               WHERE cs.biomolecule_id = $1
                     -- NO INTRAMOLECULAR OR SECONDARY CONTACTS
                     AND cs.is_same_entity = false
                     AND rbgn.chain_id != rend.chain_id
            GROUP BY cs.biomolecule_id, rbgn.chain_id, rend.chain_id
            ORDER BY 1,2,3;
            ' USING biomol_id;

            RAISE NOTICE 'inserted protein-protein interfaces for biomolecule %', biomol_id;
        END LOOP;
END$$;

    -- UPDATE THE PTREE PATHS OF THE INTERFACES
    UPDATE credo.interfaces i
       SET path = sq.path
      FROM (
            SELECT i.interface_id,
                   (
                    s.pdb || '/' ||
                    b.assembly_serial || '/' || 'I:' ||
                    cbgn.pdb_chain_id || '-' ||
                    cend.pdb_chain_id
                   ) ::ptree AS path
              FROM credo.interfaces i
              JOIN credo.biomolecules b USING(biomolecule_id)
              JOIN credo.structures s USING(structure_id)
              JOIN credo.chains cbgn ON i.chain_bgn_id = cbgn.chain_id
              JOIN credo.chains cend ON i.chain_end_id = cend.chain_id
           ) sq
     WHERE sq.interface_id = i.interface_id;

-- UPDATE INTERFACES THAT ARE ONLY IN QUATERNARY ASSEMBLIES
UPDATE credo.interfaces i
   SET is_quaternary = true
  FROM credo.chains cbgn, credo.chains cend
 WHERE cbgn.chain_id = i.chain_bgn_id AND cend.chain_id = i.chain_end_id
       AND (cbgn.is_at_identity = false OR cend.is_at_identity = false);

-- RESIDUES IN PROTEIN-PROTEIN INTERFACES
  INSERT INTO credo.interface_residues
  SELECT DISTINCT i.interface_id, rbgn.residue_id, rend.residue_id
    FROM credo.residue_interaction_pairs rip
    JOIN credo.residues rbgn ON rbgn.residue_id = rip.residue_bgn_id
    JOIN credo.residues rend ON rend.residue_id = rip.residue_end_id
    JOIN credo.interfaces i ON i.chain_bgn_id = rbgn.chain_id AND i.chain_end_id = rend.chain_id
   WHERE i.interface_id > (SELECT COALESCE(MAX(interface_id),0) FROM credo.interface_residues)
         AND rbgn.entity_type_bm & 32 != 0 AND rend.entity_type_bm & 32 != 0
         -- NOT INTRAMOLECULAR RESIDUE INTERACTIONS
         AND rbgn.chain_id != rend.chain_id
ORDER BY 1,2,3;

-- SET A FLAG FOR INTERFACES THAT CONTAIN AT LEAST ONE MODIFIED RESIDUE
UPDATE credo.interfaces i
SET    has_mod_res = true
FROM   credo.interface_residues ir, credo.peptides pbgn, credo.peptides pend
WHERE  ir.interface_id = i.interface_id
       AND pbgn.residue_id = ir.residue_bgn_id AND  pend.residue_id = ir.residue_end_id
       AND (pbgn.is_modified = true OR pend.is_modified = true);

 -- SET A FLAG FOR INTERFACES THAT CONTAIN AT LEAST ONE RESIDUE WITH MISSING ATOMS
UPDATE credo.interfaces i
SET    has_missing_atoms = true
FROM   credo.interface_residues ir, credo.peptides pbgn, credo.peptides pend
WHERE  ir.interface_id = i.interface_id
       AND pbgn.residue_id = ir.residue_bgn_id AND  pend.residue_id = ir.residue_end_id
       AND (pbgn.is_incomplete = true OR pend.is_incomplete = true);

-- POPULATE TABLE OF PROTEIN-OLIGONUCLEOTIDE GROOVES
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR biomol_id IN   SELECT biomolecule_id
                             FROM credo.biomolecules
                            WHERE biomolecule_id > (SELECT COALESCE(max(biomolecule_id),0) FROM credo.grooves)
                         ORDER BY 1
        LOOP
            -- INSERT CONTACTS
            EXECUTE
            '
               INSERT INTO credo.grooves(biomolecule_id, chain_prot_id, chain_nuc_id, num_res_prot, num_res_nuc, nucleic_acid_type)
                      -- FIRST CHAIN IS POLYPEPTIDE, SECOND IS NUCLEIC ACID
               SELECT p.biomolecule_id,
                      p.chain_id AS chain_prot_id, n.chain_id AS chain_nuc_id,
                      COUNT(DISTINCT p.residue_id) as num_res_prot,
                      COUNT(DISTINCT n.residue_id) as num_res_nuc,
                      o.nucleic_acid_type
                 FROM credo.residue_interaction_pairs rip
                 JOIN credo.peptides p ON p.residue_id = rip.residue_bgn_id
                 JOIN credo.nucleotides n ON n.residue_id = rip.residue_end_id
                 JOIN credo.oligonucleotides o ON o.chain_id = n.chain_id
                WHERE p.biomolecule_id = $1
             GROUP BY p.biomolecule_id, p.chain_id, n.chain_id, o.nucleic_acid_type
            UNION ALL
               SELECT p.biomolecule_id,
                      p.chain_id AS chain_prot_id, n.chain_id AS chain_nuc_id,
                      COUNT(DISTINCT p.residue_id) as num_res_prot,
                      COUNT(DISTINCT n.residue_id) as num_res_nuc,
                      o.nucleic_acid_type
                 FROM credo.residue_interaction_pairs rip
                 JOIN credo.peptides p ON p.residue_id = rip.residue_end_id
                 JOIN credo.nucleotides n ON n.residue_id = rip.residue_bgn_id
                 JOIN credo.oligonucleotides o ON o.chain_id = n.chain_id
                WHERE p.biomolecule_id = $1
             GROUP BY p.biomolecule_id, p.chain_id, n.chain_id, o.nucleic_acid_type
             ORDER BY 1,2,3;
            ' USING biomol_id;

            RAISE NOTICE 'inserted protein-oligonucleotide grooves for biomolecule %', biomol_id;
        END LOOP;
END$$;

-- UPDATE THE PTREE PATHS OF THE GROOVES
UPDATE credo.grooves g
   SET path = sq.path
  FROM (
        SELECT g.groove_id,
               (
                s.pdb || '/' ||
                b.assembly_serial || '/' || 'G:' ||
                cbgn.pdb_chain_id || '-' ||
                cend.pdb_chain_id
               ) ::ptree AS path
          FROM credo.grooves g
          JOIN credo.biomolecules b USING(biomolecule_id)
          JOIN credo.structures s USING(structure_id)
          JOIN credo.chains cbgn ON g.chain_prot_id = cbgn.chain_id
          JOIN credo.chains cend ON g.chain_nuc_id = cend.chain_id
       ) sq
 WHERE sq.groove_id = g.groove_id;

-- RESIDUES IN PROTEIN-DNA/RNA GROOVES
  INSERT INTO credo.groove_residues
  SELECT g.groove_id, p.residue_id, n.residue_id
    FROM credo.grooves g
    JOIN credo.peptides p ON p.chain_id = g.chain_prot_id
    JOIN credo.nucleotides n ON n.chain_id = g.chain_nuc_id
    JOIN credo.residue_interaction_pairs rip
         ON rip.residue_bgn_id = p.residue_id AND rip.residue_end_id = n.residue_id
   WHERE g.groove_id > (SELECT MAX(groove_id) FROM credo.groove_residues)
   UNION
  SELECT g.groove_id, p.residue_id, n.residue_id
    FROM credo.grooves g
    JOIN credo.peptides p ON p.chain_id = g.chain_prot_id
    JOIN credo.nucleotides n ON n.chain_id = g.chain_nuc_id
    JOIN credo.residue_interaction_pairs rip
         ON rip.residue_end_id = p.residue_id AND rip.residue_bgn_id = n.residue_id
   WHERE g.groove_id > (SELECT COALESCE(MAX(groove_id),0) FROM credo.groove_residues)
ORDER BY 1,2,3;

-- SET A FLAG FOR GROOVES THAT CONTAIN AT LEAST ONE RESIDUE WITH MISSING ATOMS
UPDATE credo.grooves g
SET    has_missing_atoms = true
FROM   credo.groove_residues gr, credo.peptides p, credo.nucleotides n
WHERE  gr.groove_id = g.groove_id
       AND p.residue_id = gr.residue_prot_id
       AND n.residue_id = gr.residue_nuc_id
       AND (p.is_incomplete = true OR n.is_incomplete = true);

-- UPDATE QUATERNARY GROOVES
 UPDATE credo.grooves g
    SET is_quaternary = true
   FROM credo.chains prot, credo.chains nuc
  WHERE prot.chain_id = g.chain_prot_id AND nuc.chain_id = g.chain_nuc_id
        AND (prot.is_at_identity = false OR nuc.is_at_identity = false);

-- BINDING SITES
  INSERT INTO credo.binding_sites(ligand_id, cath_dmns, scop_pxs, has_non_std_res,
                                  has_mod_res, has_mut_res)
  SELECT ligand_id,
         ARRAY(SELECT * FROM UNNEST(array_agg(distinct cath)) as u(e) where e IS NOT NULL) as cath_dmns,
         ARRAY(SELECT * FROM UNNEST(array_agg(distinct px)) as u(e) where e IS NOT NULL) as scop_pxs,
         bool_or(is_non_std) as has_non_std_res,
         bool_or(is_modified) as has_mod_res,
         bool_or(is_mutated) as has_mut_res
    FROM credo.binding_site_residues bsr
    JOIN credo.peptides p ON p.residue_id = bsr.residue_id
GROUP BY ligand_id
ORDER BY 1;

-- UPDATE BINDING SITE CATH SUPERFAMILY LABELS
UPDATE credo.binding_sites bs
   SET hom_superfam = sq.hom_superfam, hom_superfam_label = NULLIF(l.label, 'NULL')
  FROM (
        SELECT ligand_id, subltree(d.node, 0, 4) as hom_superfam
          FROM (
                SELECT ligand_id, unnest(cath_dmns) as dmn
                  FROM credo.binding_sites bs
               ) sq
          JOIN cath.domains d ON d.dmn = sq.dmn
       ) sq,
       cath.labels l
 WHERE sq.ligand_id = bs.ligand_id AND sq.hom_superfam = l.node;

-- SET A FLAG FOR BINDING SITES THAT HAVE MISSING ATOMS
UPDATE credo.binding_sites bs
   SET has_missing_atoms = true
  FROM credo.binding_site_residues bsr
  JOIN credo.peptides p USING(residue_id)
 WHERE bs.ligand_id = bsr.ligand_id AND p.is_incomplete = true;

-- SET A FLAG FOR BINDING SITES THAT HAVE MAPPED VARIATIONS
UPDATE credo.binding_sites bs
   SET has_mapped_var = true
  FROM credo.binding_site_residues bsr
  JOIN credo.peptides p USING(residue_id)
  JOIN variations.variation_to_pdb USING(res_map_id)
 WHERE bs.ligand_id = bsr.ligand_id;

-- set a flag for binding sites that are part of protein kinases
UPDATE credo.binding_sites bs
   SET is_kinase = true
  FROM credo.binding_site_residues bsr, credo.peptides p, credo.polypeptides pp
 WHERE bs.ligand_id = bsr.ligand_id
       AND bsr.residue_id = p.residue_id
       AND p.chain_id = pp.chain_id
       AND pp.is_kinase = true;