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
