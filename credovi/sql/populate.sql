--------------------------------------------------------------------------------
-- Start with the tables that query the raw_ tables. That makes it easier to  --
-- update CREDO because these tables will only contain the new data.          --
--------------------------------------------------------------------------------



-- INSERT STRUCTURES
  INSERT INTO credo.structures(pdb)
  SELECT DISTINCT pdb
    FROM credo.raw_atoms r
ORDER BY pdb;

-- INSERT BIOMOLECULES
  INSERT INTO credo.biomolecules(structure_id, path, assembly_serial)
  SELECT DISTINCT s.structure_id, 
                  (s.pdb || '/' || assembly_serial)::ptree as path,
                  assembly_serial
    FROM credo.structures s
    JOIN credo.raw_atoms rw on rw.pdb = s.pdb
ORDER BY structure_id, assembly_serial;

-- INSERT CHAINS
  INSERT INTO credo.chains(biomolecule_id, pdb_chain_id, pdb_chain_asu_id, path, is_at_identity)
  SELECT DISTINCT b.biomolecule_id, rw.pdb_chain_id, rw.pdb_chain_asu_id, 
                  b.path || rw.pdb_chain_id as path,
                  rw.pdb_chain_id = rw.pdb_chain_asu_id as is_at_identity
    FROM credo.raw_atoms rw 
    JOIN credo.structures s ON s.pdb = rw.pdb
    JOIN credo.biomolecules b
         ON s.structure_id = b.structure_id
         AND b.assembly_serial = rw.assembly_serial
ORDER BY b.biomolecule_id, rw.pdb_chain_id;

-- INSERT RESIDUES THROUGH BIOMOLECULE LOOP
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES WITHOUT RESIDUES
        FOR biomol_id IN   SELECT biomolecule_id 
                             FROM credo.biomolecules 
                            WHERE biomolecule_id > (SELECT max(biomolecule_id) FROM credo.residues)
                         ORDER BY 1
        LOOP
            -- EXECUTE INSERT
            EXECUTE 
            '
              INSERT INTO credo.residues(biomolecule_id, chain_id, path, res_name, res_num, 
                                         ins_code, entity_type_bm)
              SELECT DISTINCT b.biomolecule_id, c.chain_id, 
                              -- LABEL PATH
                              c.path || (rw.res_name || ''`'' || rw.res_num || TRIM(rw.ins_code))::ptree,
                              rw.res_name, rw.res_num, rw.ins_code, rw.entity_type_bm
                FROM credo.structures s
                JOIN credo.biomolecules b ON b.structure_id = s.structure_id
                JOIN credo.chains c on c.biomolecule_id = b.biomolecule_id
                JOIN credo.raw_atoms rw
                     ON s.pdb = rw.pdb
                     AND b.assembly_serial = rw.assembly_serial
                     AND c.pdb_chain_id = rw.pdb_chain_id
               WHERE b.biomolecule_id = $1
            ORDER BY 1,2,4,5
            ' USING biomol_id;

            RAISE NOTICE 'inserted residues for biomolecule %', biomol_id;
        END LOOP;
END$$;

-- INSERT ATOMS
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES WITHOUT ATOMS
        FOR biomol_id IN   SELECT biomolecule_id 
                             FROM credo.biomolecules 
                            WHERE biomolecule_id > (SELECT max(biomolecule_id) FROM credo.atoms)
                         ORDER BY 1
        LOOP
            -- INSERT ATOMS
            EXECUTE 
            '
              INSERT INTO credo.atoms(biomolecule_id, residue_id, atom_serial, group_pdb, 
                                      atom_name, alt_loc, coords, occupancy, b_factor, 
                                      element, hyb, tripos_atom_type, is_donor, is_acceptor, 
                                      is_aromatic, is_weak_acceptor, is_weak_donor, is_hydrophobe, 
                                      is_metal, is_pos_ionisable, is_neg_ionisable, is_xbond_donor, 
                                      is_xbond_acceptor, is_carbonyl_oxygen, is_carbonyl_carbon)
              SELECT b.biomolecule_id, r.residue_id, atom_serial, group_pdb, atom_name, alt_loc, 
                     coords, occupancy, b_factor, element, hyb, tripos_atom_type, is_donor, 
                     is_acceptor, is_aromatic, is_weak_acceptor, is_weak_donor, is_hydrophobe, 
                     is_metal, is_pos_ionisable, is_neg_ionisable, is_xbond_donor, is_xbond_acceptor, 
                     is_carbonyl_oxygen, is_carbonyl_carbon
                FROM credo.raw_atoms rw
                JOIN credo.structures s ON rw.pdb = s.pdb
                JOIN credo.biomolecules b 
                     ON b.structure_id = s.structure_id 
                     AND b.assembly_serial = rw.assembly_serial
                JOIN credo.chains c
                     ON c.biomolecule_id = b.biomolecule_id 
                     AND c.pdb_chain_id = rw.pdb_chain_id
                JOIN credo.residues r 
                     ON r.chain_id = c.chain_id
                     AND r.res_num = rw.res_num
                     AND r.ins_code = rw.ins_code
                     AND r.res_name = rw.res_name
               WHERE b.biomolecule_id = $1
            ORDER BY r.residue_id, atom_serial
            ' USING biomol_id;

            RAISE NOTICE 'inserted atoms for biomolecule %', biomol_id;
        END LOOP;
END$$;

-- INSERT CONTACTS THROUGH AN ANONYMOUS CODE BLOCK FOR EACH BIOMOLECULE
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES WITHOUT CONTACTS
        FOR biomol_id IN   SELECT biomolecule_id 
                             FROM credo.biomolecules 
                            WHERE biomolecule_id > (SELECT max(biomolecule_id) FROM credo.contacts)
                         ORDER BY 1
        LOOP
            -- INSERT CONTACTS
            EXECUTE 
            '
              INSERT INTO credo.contacts(biomolecule_id, atom_bgn_id, atom_end_id, distance, 
                                         structural_interaction_type_bm, is_same_entity, is_clash, is_covalent, 
                                         is_vdw_clash, is_vdw, is_proximal, is_hbond, is_weak_hbond, 
                                         is_xbond, is_ionic, is_metal_complex, is_aromatic, is_hydrophobic, is_carbonyl)
              SELECT a1.biomolecule_id, a1.atom_id, a2.atom_id, rw.distance, rw.structural_interaction_type_bm, 
                     rw.is_same_entity, is_clash, is_covalent, is_vdw_clash, is_vdw, is_proximal, is_hbond, 
                     is_weak_hbond, 
                     is_xbond, is_ionic, is_metal_complex, rw.is_aromatic, is_hydrophobic, is_carbonyl
                FROM credo.raw_contacts rw
                JOIN credo.structures s ON s.pdb = rw.pdb
                JOIN credo.biomolecules b
                     ON s.structure_id = b.structure_id
                     AND b.assembly_serial = rw.assembly_serial
                JOIN credo.atoms a1
                     ON a1.biomolecule_id = b.biomolecule_id
                     AND a1.atom_serial = rw.atom_bgn_serial
                JOIN credo.atoms a2
                     ON a2.biomolecule_id = b.biomolecule_id
                     AND a2.atom_serial = rw.atom_end_serial
               WHERE a1.biomolecule_id = $1           
            ORDER BY biomolecule_id, a1.atom_id, a2.atom_id
            ' USING biomol_id;

            RAISE NOTICE 'inserted contacts for biomolecule %', biomol_id;
        END LOOP;
END$$;

-- INSERT LIGANDS
  INSERT INTO credo.ligands(biomolecule_id, path, entity_serial, pdb_chain_id, res_num, 
                            ligand_name, num_hvy_atoms)
  SELECT DISTINCT b.biomolecule_id,
                  -- PTREE LABEL
                  b.path || 
                  (
                   rw.pdb_chain_id || '/' 
                   || rw.ligand_name ||
                   COALESCE('`' || rw.res_num, '')
                  )::ptree,
                  entity_serial, rw.pdb_chain_id, rw.res_num, 
                  rw.ligand_name, rw.num_hvy_atoms
    FROM credo.structures s
    JOIN credo.biomolecules b ON b.structure_id = s.structure_id
    JOIN credo.raw_ligands rw
         ON rw.pdb = s.pdb
         AND rw.assembly_serial = b.assembly_serial
ORDER BY 1,3,4;

-- INSERT LIGAND_COMPONENTS // SLOW SHOULD BE REWRITTEN
  INSERT INTO credo.ligand_components(ligand_id, residue_id)
  SELECT DISTINCT l.ligand_id, r.residue_id
    FROM credo.raw_atoms rw
    JOIN credo.structures s ON s.pdb = rw.pdb
    JOIN credo.biomolecules b
         ON b.structure_id = s.structure_id
         AND b.assembly_serial = rw.assembly_serial
    JOIN credo.chains c
         ON c.biomolecule_id = b.biomolecule_id
         AND c.pdb_chain_id = rw.pdb_chain_id
    JOIN credo.residues r
         ON c.chain_id = r.chain_id
         AND r.res_num = rw.res_num
         AND r.ins_code = rw.ins_code
         AND r.res_name = rw.res_name
    JOIN credo.ligands l
         ON l.biomolecule_id = b.biomolecule_id
         AND l.entity_serial = rw.entity_serial
ORDER BY l.ligand_id, r.residue_id;

-- INSERT AROMATIC RINGS
  INSERT INTO credo.aromatic_rings(biomolecule_id, residue_id, ring_serial, centroid, size)
  SELECT b.biomolecule_id, a.residue_id, ring_serial,
         AVG(coords) as centroid,
         COUNT(rw.atom_serial) as size
    FROM credo.raw_aromaticrings rw
    JOIN credo.structures s on s.pdb = rw.pdb
    JOIN credo.biomolecules b
         ON b.structure_id = s.structure_id
         AND rw.assembly_serial = b.assembly_serial
    JOIN credo.atoms a
         ON a.biomolecule_id = b.biomolecule_id AND a.atom_serial = rw.atom_serial
GROUP BY b.biomolecule_id, residue_id, ring_serial
ORDER BY b.biomolecule_id, residue_id, ring_serial;

        -- UPDATE THE AROMATIC RING SERIAL NUMBER OF THE RESIDUE
        -- REQUIRED FOR RING-INTERACTIONS
        UPDATE credo.aromatic_rings ar
           SET ring_number = rs.ring_number
          FROM (
               SELECT ar.aromatic_ring_id, ar.residue_id, ar.biomolecule_id,
                      ROW_NUMBER() OVER (PARTITION BY residue_id ORDER BY residue_id) AS ring_number
                 FROM credo.aromatic_rings ar
               ) rs
         WHERE ar.aromatic_ring_id = rs.aromatic_ring_id
               AND ar.ring_number IS NULL;

        -- UPDATE PATHS OF AROMATIC RINGS
        UPDATE credo.aromatic_rings ar
           SET path = r.path || ('AR:' || ring_number)::ptree
          FROM credo.residues r 
         WHERE r.residue_id = ar.residue_id;

-- AROMATIC RING ATOMS
  INSERT INTO credo.aromatic_ring_atoms(aromatic_ring_id, atom_id)
  SELECT DISTINCT aromatic_ring_id, atom_id
    FROM credo.aromatic_rings ar
    JOIN credo.atoms a on a.residue_id = ar.residue_id
    JOIN credo.biomolecules b ON b.biomolecule_id = a.biomolecule_id
    JOIN credo.structures s on s.structure_id = b.structure_id
    JOIN credo.raw_aromaticrings rw
         ON rw.pdb = s.pdb
         AND rw.assembly_serial = b.assembly_serial
         AND rw.ring_serial = ar.ring_serial
         AND rw.atom_serial = a.atom_serial
ORDER BY aromatic_ring_id, atom_id;

        -- SET THE NORMAL FOR EACH AROMATIC RING
        -- REQUIRED TO CALCULATE RING INTERACTIONS
        UPDATE credo.aromatic_rings ar
           SET normal = rs.normal
          FROM (
                 SELECT ar.residue_id, ring_serial, three_point_normal(array_agg(coords)) as normal
                   FROM credo.aromatic_rings ar
                   JOIN credo.aromatic_ring_atoms ara ON ar.aromatic_ring_id = ara.aromatic_ring_id
                   JOIN credo.atoms a ON a.atom_id = ara.atom_id
                  WHERE ar.normal IS NULL
               GROUP BY ar.residue_id, ring_serial
               ORDER BY ar.residue_id, ring_serial
               ) rs
         WHERE rs.residue_id = ar.residue_id AND rs.ring_serial = ar.ring_serial;



--------------------------------------------------------------------------------
-- Continue with the tables that do not depend on the raw_ tables. The update --
-- process is slighty more complicated because either the tables have to be   --
-- truncated or LEFT JOINs used to exclude already existing values.           --
--------------------------------------------------------------------------------



-- INSERT NEW HETATMS
   INSERT INTO credo.hetatms(atom_id, ligand_component_id, ligand_id)
   SELECT a.atom_id, lc.ligand_component_id, lc.ligand_id
     FROM credo.atoms a
     JOIN credo.residues r ON r.residue_id = a.residue_id
     JOIN credo.ligand_components lc ON lc.residue_id = r.residue_id
LEFT JOIN credo.hetatms h ON h.atom_id = a.atom_id
    WHERE h.atom_id IS NULL;

-- CREATE LIGAND TO BINDING-SITE RESIDUE MAPPING
-- CAN BE SAFELY TRUNCATED DUE TO LACK OF SERIAL COLUMN
TRUNCATE TABLE credo.binding_sites;
  INSERT INTO credo.binding_sites
  SELECT h.ligand_id, p.residue_id
    FROM credo.hetatms h
    JOIN credo.contacts cs on h.atom_id = cs.atom_bgn_id
    JOIN credo.atoms a on cs.atom_end_id = a.atom_id
    JOIN credo.peptides p ON p.residue_id = a.residue_id
   WHERE a.biomolecule_id = cs.biomolecule_id 
         AND cs.is_same_entity = false
   UNION
  SELECT h.ligand_id, p.residue_id
    FROM credo.hetatms h
    JOIN credo.contacts cs on h.atom_id = cs.atom_end_id
    JOIN credo.atoms a on cs.atom_bgn_id = a.atom_id
    JOIN credo.peptides p ON p.residue_id = a.residue_id
   WHERE a.biomolecule_id = cs.biomolecule_id  
         AND cs.is_same_entity = false
ORDER BY 1,2;

-- CALCULATE RING INTERACTIONS
  INSERT INTO credo.ring_interactions(biomolecule_id,
                                      aromatic_ring_bgn_id, aromatic_ring_end_id, 
                                      distance, dihedral, theta, iota)
  SELECT ar1.biomolecule_id,
         ar1.aromatic_ring_id, ar2.aromatic_ring_id,
         ar1.centroid -> ar2.centroid AS distance,
         SIGNED_DEGREES(ar1.normal @ ar2.normal) AS dihedral,
         SIGNED_DEGREES(ar1.normal @ (ar1.centroid - ar2.centroid)) AS theta,
         SIGNED_DEGREES(ar2.normal @ (ar2.centroid - ar1.centroid)) AS iota
    FROM credo.aromatic_rings ar1, credo.aromatic_rings ar2
   WHERE ar1.aromatic_ring_id < ar2.aromatic_ring_id
         -- CONSIDER ALL PAIRS INSIDE THE SAME STRUCTURE
         AND ar1.biomolecule_id = ar2.biomolecule_id
         -- MAXIMUM DISTANCE BETWEEN CENTROIDS
         AND ar1.centroid -> ar2.centroid < 6.0
         -- RINGS HAVE TO BE FROM DIFFERENT RESIDUES
         AND ar2.residue_id != ar1.residue_id
         --
         AND ar1.biomolecule_id > (SELECT MAX(biomolecule_id) FROM credo.ring_interactions)
ORDER BY 1,2;

-- INSERT ATOM-AROMATIC RING INTERACTIONS
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES THAT DO NOT HAVE ANY RING INTERACTIONS YET
        FOR biomol_id IN      SELECT biomolecule_id
                                FROM (
                                      SELECT biomolecule_id FROM credo.biomolecules
                                      EXCEPT 
                                      SELECT biomolecule_id 
                                        FROM credo.aromatic_rings 
                                        JOIN credo.atom_ring_interactions USING(aromatic_ring_id)
                                     ) sq
                            ORDER BY 1
        LOOP
            EXECUTE 
                '
                INSERT      INTO credo.atom_ring_interactions(aromatic_ring_id, atom_id, distance, theta, interaction_type)
                SELECT      DISTINCT ar.aromatic_ring_id, a.atom_id,
                            ar.centroid -> a.coords as distance,
                            SIGNED_DEGREES(ar.normal @ (ar.centroid - a.coords)) AS theta,
                            CASE
                                WHEN a.element = ''C'' AND a.is_weak_donor = true THEN ''CARBONPI''
                                WHEN a.is_pos_ionisable THEN ''CATIONPI''
                                WHEN a.is_donor THEN ''DONORPI''
                                WHEN a.is_xbond_donor THEN ''HALOGENPI''
                                ELSE NULL
                            END AS interaction_type
                FROM        credo.aromatic_rings ar, credo.atoms a
                WHERE       ar.residue_id != a.residue_id 
                            AND a.biomolecule_id = ar.biomolecule_id
                            AND a.is_aromatic = false
                            AND (ar.centroid -> a.coords) <= 4.5
                            AND (SELECT DEGREES(ar.normal @+ (ar.centroid - a.coords)) <= 30.0)
                            AND a.biomolecule_id = $1
                ORDER BY    ar.aromatic_ring_id, a.atom_id;
                '
            USING biomol_id;

            RAISE NOTICE 'Atom-aromatic ring interactions updated for Biomolecule %', biomol_id;
        END LOOP;
END$$;

-- INSERT CHAIN PROTEIN FRAGMENTS AND FRAGMENT TO RESIDUE MAPPING
DO $$
    DECLARE
        cur_chain_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL CHAINS THAT ARE NOT IN THE PROTFRAGMENTS TABLE YET
        FOR cur_chain_id IN   SELECT chain_id 
                                FROM credo.chains 
                               WHERE chain_id > (SELECT max(chain_id) FROM credo.prot_fragments)
                            ORDER BY 1
        LOOP
            -- INSERT PROTEIN SECONDARY STRUCTURE FRAGMENTS
            EXECUTE 
                '
                    INSERT INTO credo.prot_fragments(chain_id, path, sstruct_serial, sstruct, 
                                                     fragment_size, fragment_seq)
                    SELECT DISTINCT 
                           c.chain_id, 
                           c.path || (''PF:'' ||  f.sstruct_serial)::ptree,
                           f.sstruct_serial, f.sstruct, f.fragment_size, f.fragment_seq
                     FROM credo.structures s
                          JOIN credo.biomolecules b ON b.structure_id = s.structure_id
                          JOIN credo.chains c on c.biomolecule_id = b.biomolecule_id
                     JOIN pdb.res_map m
                          ON m.pdb = s.pdb AND m.pdb_chain_id = c.pdb_chain_id
                     JOIN pdb.pdb_prot_fragments f
                          ON f.pdb = m.pdb
                          AND f.pdb_chain_id = m.pdb_chain_id
                          AND f.sstruct_serial = m.sstruct_serial
                    WHERE c.chain_id = $1
                 ORDER BY c.chain_id, f.sstruct_serial
                ' USING cur_chain_id;

            -- UPDATE THE MAPPING BETWEEN PROTEIN FRAGMENTS AND THEIR RESIDUES
            EXECUTE
                '
                   INSERT INTO credo.prot_fragment_residues
                   SELECT DISTINCT pf.prot_fragment_id, p.residue_id
                     FROM credo.prot_fragments pf
                     JOIN credo.chains c ON c.chain_id = pf.chain_id
                     JOIN credo.biomolecules b ON c.biomolecule_id = b.biomolecule_id
                     JOIN credo.structures s ON s.structure_id = b.structure_id
                     JOIN credo.peptides p ON p.chain_id = c.chain_id
                     JOIN pdb.res_map m
                          ON m.pdb = s.pdb
                          AND m.pdb_chain_id = c.pdb_chain_asu_id
                          AND m.pdb_res_num = p.res_num
                          AND m.pdb_ins_code = p.ins_code
                          AND pf.sstruct_serial = m.sstruct_serial
                    WHERE c.chain_id = $1
                 ORDER BY 1,2;
                ' USING cur_chain_id;

            RAISE NOTICE 'inserted protein fragments and update residue mapping for chain %', cur_chain_id;
        END LOOP;
END$$;

-- UPDATE PROTEIN-FRAGMENT COMPLETENESS
UPDATE credo.prot_fragments pf
   SET completeness = rs.residues / pf.fragment_size::real
  FROM (
          SELECT pfr.prot_fragment_id, count(residue_id) as residues
            FROM credo.prot_fragment_residues pfr
        GROUP BY pfr.prot_fragment_id
       ) rs
 WHERE rs.prot_fragment_id = pf.prot_fragment_id;


-- RESIDUE INTERACTION PAIRS
-- THIS TABLE IS MOSTLY USED TO UPDATE THE INTERFACE/GROOVE RESIDUES TABLES
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR biomol_id IN   SELECT biomolecule_id 
                             FROM credo.biomolecules 
                            WHERE biomolecule_id > (SELECT max(biomolecule_id) FROM credo.residue_interaction_pairs)
                         ORDER BY 1
        LOOP
            -- INSERT CONTACTS
            EXECUTE 
                '
                  INSERT INTO credo.residue_interaction_pairs
                  SELECT DISTINCT cs.biomolecule_id, abgn.residue_id, aend.residue_id, 
                         cs.structural_interaction_type_bm
                    FROM credo.contacts cs
                    JOIN credo.atoms abgn 
                         ON abgn.atom_id = cs.atom_bgn_id AND cs.biomolecule_id = abgn.biomolecule_id
                    JOIN credo.atoms aend 
                         ON aend.atom_id = cs.atom_end_id AND cs.biomolecule_id = aend.biomolecule_id
                   WHERE cs.biomolecule_id = $1
                ORDER BY 1,2
                ' USING biomol_id;

            RAISE NOTICE 'inserted residue interaction pairs for biomolecule %', biomol_id;
        END LOOP;
END$$;



------------------------------------
-- UPDATE COLUMNS                 --
-- ONLY UPDATES FROM CREDO TABLES --
------------------------------------



-- UPDATE NUMBER OF BIOMOLECULES PER STRUCTURE
UPDATE      credo.structures s
SET         num_biomolecules = sq.biomolecules
FROM        (
            SELECT      structure_id, COUNT(biomolecule_id) as biomolecules
            FROM        credo.biomolecules b
            GROUP BY    structure_id
            ) sq
WHERE       s.structure_id = sq.structure_id;

-- UPDATE BIOMOLECULE ASSEMBLY TYPE
UPDATE      credo.biomolecules b
SET         assembly_type = rs.assembly_type
FROM        (
            SELECT  sq.biomolecule_id,
                    CASE
                        WHEN num_chains = 1 THEN 'monomer'
                        WHEN num_chains = 2 THEN 'dimer'
                        WHEN num_chains = 3 THEN 'trimer'
                        WHEN num_chains = 4 THEN 'tetramer'
                        WHEN num_chains = 5 THEN 'pentamer'
                        WHEN num_chains = 6 THEN 'hexamer'
                        WHEN num_chains = 7 THEN 'heptamer'
                        WHEN num_chains = 8 THEN 'octamer'
                        WHEN num_chains = 9 THEN 'nonamer'
                        WHEN num_chains = 10 THEN 'decamer'
                        WHEN num_chains = 11 THEN 'undecamer'
                        WHEN num_chains = 12 THEN 'dodecamer'
                        WHEN num_chains = 13 THEN 'tridecamer'
                        WHEN num_chains = 14 THEN 'tetradecamer'
                        WHEN num_chains = 15 THEN 'pentadecamer'
                        WHEN num_chains = 16 THEN 'hexadecamer'
                        WHEN num_chains = 17 THEN 'heptadecamer'
                        WHEN num_chains = 18 THEN 'octadecamer'
                        WHEN num_chains = 19 THEN 'nonadecamer'
                        WHEN num_chains = 20 THEN 'eicosamer'
                        ELSE num_chains || '-mer'
                    END as assembly_type
            FROM    (
                    SELECT      biomolecule_id,
                                COUNT(DISTINCT chain_id) as num_chains
                    FROM        credo.chains c
                    GROUP BY    biomolecule_id
                    ) sq
            ) rs
WHERE       b.biomolecule_id = rs.biomolecule_id;

-- UPDATE THE NUMBER OF LIGANDS  FOR EACH BIOMOLECULE
UPDATE      credo.biomolecules b
SET         num_ligands = rs.num_ligands
FROM        (
            SELECT      biomolecule_id,
                        COUNT(DISTINCT ligand_id) as num_ligands
            FROM        credo.ligands l
            GROUP BY    biomolecule_id
            ) rs
WHERE       b.biomolecule_id = rs.biomolecule_id;

-- UPDATE THE NUMBER OF CHAINS AND ATOMS FOR EACH BIOMOLECULE
UPDATE      credo.biomolecules b
SET         num_chains = rs.num_chains, num_atoms = rs.num_atoms
FROM        (
            SELECT      r.biomolecule_id,
                        COUNT(DISTINCT chain_id) as num_chains,
                        COUNT(DISTINCT atom_id) as num_atoms
            FROM        credo.residues r
            JOIN        credo.atoms a ON a.residue_id = r.residue_id
            WHERE       a.biomolecule_id = r.biomolecule_id
            GROUP BY    r.biomolecule_id
            ) rs
WHERE       b.biomolecule_id = rs.biomolecule_id;

-- SET A FLAG FOR CHAINS THAT HAVE DISORDERED REGIONS
UPDATE credo.chains c
   SET has_disordered_regions = true
  FROM credo.biomolecules b, credo.structures s, pdb.disordered_regions dr
 WHERE c.biomolecule_id = b.biomolecule_id 
       AND b.structure_id = s.structure_id 
       AND dr.pdb = s.pdb AND dr.pdb_chain_id = c.pdb_chain_asu_id;

-- SET A FLAG FOR DISORDERED RESIDUES
UPDATE      credo.residues r
SET         is_disordered = TRUE
FROM        credo.atoms a
WHERE       a.residue_id = r.residue_id AND a.alt_loc != ' ' AND a.biomolecule_id = r.biomolecule_id;

-- SET A FLAG FOR DISORDERED LIGANDS
UPDATE      credo.ligands l
SET         is_disordered = true
FROM        credo.hetatms h, credo.atoms a
WHERE       l.ligand_id = h.ligand_id AND h.atom_id = a.atom_id 
            AND a.alt_loc != ' ';

-- SET A FLAG IF THE LIGAND IS FROM ASU OR NOT
UPDATE      credo.ligands l
SET         is_at_identity = true
FROM        credo.biomolecules b, credo.chains c
WHERE       l.biomolecule_id = b.biomolecule_id 
            AND b.biomolecule_id = c.biomolecule_id 
            AND l.pdb_chain_id = c.pdb_chain_id
            AND c.is_at_identity = true;

-- UPDATE N- AND C-TERMINUS INFORMATION OF PROTEIN-FRAGMENTS
UPDATE      credo.prot_fragments pf
SET         prot_fragment_nterm_id = sq.prot_fragment_nterm_id,
            prot_fragment_cterm_id = sq.prot_fragment_cterm_id
FROM        (
            SELECT  prot_fragment_id,
                    LAG(prot_fragment_id) OVER (PARTITION BY chain_id ORDER BY chain_id) AS prot_fragment_nterm_id,
                    LEAD(prot_fragment_id) OVER (PARTITION BY chain_id ORDER BY chain_id) AS prot_fragment_cterm_id
            FROM    credo.prot_fragments pf
            ) sq
WHERE       pf.prot_fragment_id = sq.prot_fragment_id;

-- UPDATE HETEROAROMATIC FLAG
UPDATE      credo.aromatic_rings ar
SET         is_hetero_aromatic = true
FROM        (
            SELECT      DISTINCT ar.aromatic_ring_id
            FROM        credo.aromatic_rings ar
            JOIN        credo.aromatic_ring_atoms ara
                        ON ar.aromatic_ring_id = ara.aromatic_ring_id
            JOIN        credo.atoms a ON a.atom_id = ara.atom_id
            WHERE       a.element != 'C'
            ) rs
WHERE       rs.aromatic_ring_id = ar.aromatic_ring_id;

-- SET THE RING INTERACTION TYPE FOR ALL RING INTERACTIONS
UPDATE      credo.ring_interactions ri
SET         interaction_type = rs1.interaction_type
FROM        (
            SELECT      ring_interaction_id,
                        CASE
                            WHEN dihedral <= 30.0 AND theta <= 30.0 THEN 'FF'
                            WHEN dihedral <= 30.0 AND theta <= 60.0 THEN 'OF'
                            WHEN dihedral <= 30.0 AND theta <= 90.0 THEN 'EE'

                            WHEN dihedral > 30.0 AND dihedral <= 60.0 AND theta <= 30.0 THEN 'FT'
                            WHEN dihedral > 30.0 AND dihedral <= 60.0 AND theta <= 60.0 THEN 'OT'
                            WHEN dihedral > 30.0 AND dihedral <= 60.0 AND theta <= 90.0 THEN 'ET'

                            WHEN dihedral > 60.0 AND dihedral <= 90.0 AND theta <= 30.0 THEN 'FE'
                            WHEN dihedral > 60.0 AND dihedral <= 90.0 AND theta <= 60.0 THEN 'OE'
                            WHEN dihedral > 60.0 AND dihedral <= 90.0 AND theta <= 90.0 THEN 'EF'
                        END as interaction_type
            FROM        (
                        SELECT  ring_interaction_id, ABS(dihedral) as dihedral, ABS(theta) as theta
                        FROM    credo.ring_interactions
                        ) rs2
            ) rs1
WHERE       ri.ring_interaction_id = rs1.ring_interaction_id;

-- CLOSEST ATOM DETAILS FOR RING INTERACTIONS
UPDATE      credo.ring_interactions ri
SET         closest_atom_bgn_id = sq.atom_bgn_id, closest_atom_End_id = sq.atom_end_id, closest_atom_distance = sq.closest_atom_distance
FROM        (
            WITH sq AS
            (
                SELECT  ri.ring_interaction_id, a1.atom_id AS atom_bgn_id, a2.atom_id as atom_end_id,
                        a1.coords -> a2.coords as closest_atom_distance,
                        ROW_NUMBER() OVER(PARTITION BY ri.ring_interaction_id
                                        ORDER BY ri.ring_interaction_id, a1.coords -> a2.coords) AS rank
                FROM    credo.ring_interactions ri
                JOIN    credo.aromatic_ring_atoms ara1 ON ara1.aromatic_ring_id = ri.aromatic_ring_bgn_id
                JOIN    credo.aromatic_ring_atoms ara2 ON ara2.aromatic_ring_id = ri.aromatic_ring_end_id
                JOIN    credo.atoms a1 ON a1.atom_id = ara1.atom_id
                JOIN    credo.atoms a2 ON a2.atom_id = ara2.atom_id
            )
            SELECT  sq.ring_interaction_id, sq.atom_bgn_id, sq.atom_end_id, sq.closest_atom_distance
            FROM    sq
            WHERE   rank = 1
            ) sq
WHERE       sq.ring_interaction_id = ri.ring_interaction_id;

-- UPDATE SECONDARY CONTACTS
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR biomol_id IN SELECT biomolecule_id 
                           FROM credo.biomolecules 
                          WHERE biomolecule_id > (SELECT COALESCE(max(biomolecule_id),0)
                                                    FROM credo.contacts
                                                   WHERE is_secondary = true)
        LOOP
            -- INSERT CONTACTS
            EXECUTE 
                '
                UPDATE      credo.contacts cs
                SET         is_secondary = true
                FROM        (
                            SELECT  cs.contact_id
                            FROM    credo.contacts cs
                                    -- GET THE COORDINATES OF BOTH ATOMS THAT FORM THIS CONTACT PAIR
                            JOIN    credo.atoms A 
                                    ON A.atom_id = cs.atom_bgn_id AND A.biomolecule_id = cs.biomolecule_id
                            JOIN    credo.atoms B 
                                    ON B.atom_id = cs.atom_end_id AND B.biomolecule_id = cs.biomolecule_id
                                    -- GET ALL POSSIBLE THIRD ATOMS THAT ARE IN CONTACT WITH THE FIRST ATOM (BUT NOT SAME RESIDUE!)
                            JOIN    credo.contacts cs2 
                                    ON cs2.atom_bgn_id = A.atom_id AND cs2.biomolecule_id = A.biomolecule_id
                            JOIN    credo.atoms C ON C.atom_id = cs2.atom_end_id
                            WHERE   B.atom_id != C.atom_id
                                    -- THE THIRD ATOM HAS TO BE CLOSER TO THE FIRST
                                    AND cs2.distance < cs.distance
                                    -- THE DISTANCE BETWEEN THE SECOND AND THIRD HAS TO BE SMALLER THAN THE DISTANCE BETWEEN FIRST AND SECOND
                                    AND (B.coords -> C.coords) < cs.distance
                                    -- THE TRIANGLE HEIGHT
                                    AND cos_triangle_height(B.coords -> C.coords, cs2.distance, cs.distance) <= 1.2
                                    AND A.biomolecule_id = $1
                            UNION
                            SELECT  cs.contact_id
                            FROM    credo.contacts cs
                                    -- GET THE COORDINATES OF BOTH ATOMS THAT FORM THIS CONTACT PAIR
                            JOIN    credo.atoms A 
                                    ON A.atom_id = cs.atom_end_id 
                                    AND A.biomolecule_id = cs.biomolecule_id
                            JOIN    credo.atoms B 
                                    ON B.atom_id = cs.atom_bgn_id 
                                    AND B.biomolecule_id = cs.biomolecule_id
                                    -- GET ALL POSSIBLE THIRD ATOMS THAT ARE IN CONTACT WITH THE FIRST ATOM (BUT NOT SAME RESIDUE!)
                            JOIN    credo.contacts cs2 
                                    ON cs2.atom_end_id = A.atom_id 
                                    AND cs2.biomolecule_id = A.biomolecule_id
                            JOIN    credo.atoms C ON C.atom_id = cs2.atom_bgn_id
                            WHERE   B.atom_id != C.atom_id
                                    -- THE THIRD ATOM HAS TO BE CLOSER TO THE FIRST
                                    AND cs2.distance < cs.distance
                                    -- THE DISTANCE BETWEEN THE SECOND AND THIRD HAS TO BE SMALLER THAN THE DISTANCE BETWEEN FIRST AND SECOND
                                    AND (B.coords -> C.coords) < cs.distance
                                    -- THE TRIANGLE HEIGHT
                                    AND cos_triangle_height(B.coords -> C.coords, cs2.distance, cs.distance) <= 1.2
                                    AND A.biomolecule_id = $1
                            ) rs
                WHERE       rs.contact_id = cs.contact_id AND cs.biomolecule_id = $1;
                ' USING biomol_id;

            RAISE NOTICE 'Updated secondary contacts for Biomolecule %', biomol_id;
        END LOOP;
END$$;

/*
This anonymous code block updates the Gini coefficient based on the number of atoms that are in contact
for each ligand with at least 5 heavy atoms.
*/
DO $$
    # GET ALL LIGAND IDS
    result = plpy.execute("SELECT ligand_id FROM credo.ligands WHERE num_hvy_atoms >= 7 AND gini_index_contacts IS NULL ORDER BY 1")

    for ligand_id in [row['ligand_id'] for row in result]:
        
        # GET ALL THE PRIMARY CONTACTS FOR EACH LIGAND
        statement = '''
                    SELECT  array_length(array_agg(sq.atoms),1) as num_atoms
                    FROM    (
                            SELECT  h.ligand_id, h.hetatm_id, cs.atom_end_id as atoms
                            FROM    credo.contacts cs
                            JOIN    credo.hetatms h ON cs.atom_bgn_id = h.atom_id
                            JOIN    credo.ligands l ON l.ligand_id = h.ligand_id
                            WHERE   l.ligand_id = $1 
                                    AND l.biomolecule_id = cs.biomolecule_id 
                                    AND cs.is_secondary = false
                            UNION
                            SELECT  h.ligand_id, h.hetatm_id, cs.atom_bgn_id as atoms
                            FROM    credo.contacts cs
                            JOIN    credo.hetatms h ON cs.atom_end_id = h.atom_id
                            JOIN    credo.ligands l ON l.ligand_id = h.ligand_id
                            WHERE   h.ligand_id = $1 
                                    AND l.biomolecule_id = cs.biomolecule_id 
                                    AND cs.is_secondary = false
                            ) sq
                    GROUP BY ligand_id, hetatm_id
                    ORDER BY 1
                    '''
                    
        # PREPARE STATEMENT FOR EXECUTION            
        query = plpy.prepare(statement, ['integer'])
        update = plpy.prepare("UPDATE credo.ligands SET gini_index_contacts = $2 WHERE ligand_id = $1", ['integer','float'])

        # FETCH THE CONTACTS
        result = plpy.execute(query, [ligand_id,])

        # TURN RESULTSET INTO LIST
        contacts = [row['num_atoms'] for row in result]

        # TOTAL NUMBER OF CONTACTS
        N = len(contacts)

        if not N: 
            plpy.warning("No primary contacts found for ligand %i"  % (ligand_id))
            continue

        # CALCULATE GINI INDEX
        try:
            gini = sum(contacts[i] * (N-i) for i in xrange(N)) 
            gini = (2.0 * gini) / (N * sum(contacts))
            gini = (1 + 1.0 / N) - gini
        except ZeroDivisionError:
            plpy.warning("cannot calculate Gini Index for ligand %i: ZeroDivisionError"  % (ligand_id))
            gini = 0.0

        # UPDATE TABLE
        plpy.execute(update, [ligand_id, gini])
$$ LANGUAGE plpythonu;
