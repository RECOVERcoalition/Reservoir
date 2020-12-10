SELECT md.molregno,
       md.pref_name,
       md.chembl_id,
       md.max_phase,
       md.molecule_type,
       cs.canonical_smiles
FROM molecule_dictionary md
LEFT JOIN compound_structures cs ON md.molregno=cs.molregno;
