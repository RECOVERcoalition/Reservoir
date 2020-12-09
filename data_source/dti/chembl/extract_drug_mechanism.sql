SELECT md.pref_name,
       md.max_phase,
       md.chembl_id AS chembl_molecule_id,
       cs.canonical_smiles,
       cs.standard_inchi,
       cs.standard_inchi_key,
       dm.action_type,
       dm.MECHANISM_OF_ACTION,
       dm.direct_interaction,
       td.chembl_id AS chembl_target_id,
       td.organism
FROM compound_records cr,
     molecule_dictionary md,
     compound_structures cs,
     drug_mechanism dm,
     target_dictionary td
WHERE md.molregno=dm.molregno
  AND dm.tid = td.tid
  AND md.molregno=cs.molregno
  AND cr.molregno=cs.molregno;
