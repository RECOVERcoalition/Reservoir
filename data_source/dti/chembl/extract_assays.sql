SELECT md.chembl_id AS chembl_molecule_id,
       cs.canonical_smiles,
       td.chembl_id AS chembl_target_id,
       td.organism,
       ac.pchembl_value,
       ac.standard_type,
       ass.assay_type
FROM compound_records cr,
     compound_structures cs,
     activities ac,
     assays ass,
     target_dictionary td,
     docs,
     molecule_dictionary md
WHERE md.molregno=cs.molregno
  AND cr.record_id=ac.record_id
  AND cr.molregno=cs.molregno
  AND ac.doc_id=docs.doc_id
  AND ac.assay_id=ass.assay_id
  AND ass.tid=td.tid
  AND ass.assay_type in ('A',
                         'F',
                         'B')
  AND NOT ac.pchembl_value IS NULL
  AND NOT td.organism IS NULL;
