import pandas as pd
import reservoir as rsv
from reservoir.parsers.proteins import protein_mapper
from reservoir.parsers.drugs import molecule_mapper

# load dtis
chembl_dtis = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/dti/chembl/all_assays_abridged.csv"
)
chembl_dm = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/dti/chembl/drug_mechanism.csv"
)

""" Map targets to hgnc
"""
# chembl id to uniprot
chembl_target_to_uniprot = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/dti/chembl/chembl_uniprot_mapping.tsv",
    sep="\t",
    names=["uniprot_id", "chembl_id", "protein_name", "protein_type"],
)

# map the uniprot id to a gene for dtis
chembl_dtis = chembl_dtis.merge(
    chembl_target_to_uniprot, left_on="chembl_target_id", right_on="chembl_id"
)
chembl_dtis["gene_hgnc_id"] = chembl_dtis["uniprot_id"].apply(
    protein_mapper.uniprot_to_hgnc
)
chembl_dtis = chembl_dtis.loc[
    ~pd.isnull(chembl_dtis["gene_hgnc_id"]),
    [
        "chembl_molecule_id",
        "gene_hgnc_id",
        "organism",
        "pchembl_value",
        "standard_type",
        "assay_type",
    ],
]

# map the uniprot id to a gene for dm
chembl_dm = chembl_dm.merge(
    chembl_target_to_uniprot, left_on="chembl_target_id", right_on="chembl_id"
)
chembl_dm["gene_hgnc_id"] = chembl_dm["uniprot_id"].apply(
    protein_mapper.uniprot_to_hgnc
)
chembl_dm = chembl_dm.loc[
    ~pd.isnull(chembl_dm["gene_hgnc_id"]), ["chembl_molecule_id", "gene_hgnc_id"]
]

""" Map the drug to a recover ID
"""
# get drugs from dtis and dms
chembl_drugs_from_dti = chembl_dtis[["chembl_molecule_id"]].rename(
    columns={"chembl_molecule_id": "chembl_id"}
)
chembl_drugs_from_dm = chembl_dm[["chembl_molecule_id"]].rename(
    columns={"chembl_molecule_id": "chembl_id"}
)
chembl_drugs = pd.concat(
    [chembl_drugs_from_dti, chembl_drugs_from_dm]
).drop_duplicates()

# map drugs to recover id
chembl_drugs_post_mapping = molecule_mapper.map_drugs(chembl_drugs)[
    ["chembl_id", "recover_id"]
]

# connect to dtis and dms
chembl_dtis = chembl_dtis.merge(
    chembl_drugs_post_mapping, left_on="chembl_molecule_id", right_on="chembl_id"
)[
    [
        "recover_id",
        "gene_hgnc_id",
        "organism",
        "pchembl_value",
        "standard_type",
        "assay_type",
    ]
]
chembl_dm = chembl_dm.merge(
    chembl_drugs_post_mapping, left_on="chembl_molecule_id", right_on="chembl_id"
)[["recover_id", "gene_hgnc_id"]]

# export
chembl_dtis.to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/dti/chembl_dtis.csv", index=False
)
chembl_dm.to_csv(rsv.RESERVOIR_DATA_FOLDER + "/parsed/dti/chembl_dm.csv", index=False)
