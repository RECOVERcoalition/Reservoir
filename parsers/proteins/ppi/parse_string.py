import pandas as pd
import reservoir as rsv
from reservoir.parsers.proteins import protein_mapper
import ppi_helper

# map from ensembl protein id to ensembl gene or hgnc symbol
string_protein_ids = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/ppi/string_human_v11_protein_ids.txt", sep="\t"
)
string_protein_ids["gene_hgnc_id"] = string_protein_ids["preferred_name"].apply(
    protein_mapper.hgnc_normalize
)
string_protein_ids.loc[
    pd.isnull(string_protein_ids["gene_hgnc_id"]), "gene_hgnc_id"
] = string_protein_ids.loc[
    pd.isnull(string_protein_ids["gene_hgnc_id"]), "preferred_name"
].apply(
    protein_mapper.ensembl_to_hgnc
)
string_protein_ids = string_protein_ids.loc[
    ~pd.isnull(string_protein_ids["gene_hgnc_id"])
]

# map the ids in the PPI
string = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/ppi/string_human_v11.txt", sep=" "
)
string = string.merge(
    string_protein_ids[["protein_external_id", "gene_hgnc_id"]],
    left_on="protein1",
    right_on="protein_external_id",
).rename(columns={"gene_hgnc_id": "gene_1_hgnc_id"})
string = string.merge(
    string_protein_ids[["protein_external_id", "gene_hgnc_id"]],
    left_on="protein2",
    right_on="protein_external_id",
).rename(columns={"gene_hgnc_id": "gene_2_hgnc_id"})
string = string[["gene_1_hgnc_id", "gene_2_hgnc_id", "combined_score"]]

# remove duplicate edges
string_full = ppi_helper.remove_duplicate_edges(string)
string_high_confidence = ppi_helper.remove_duplicate_edges(
    string.loc[string["combined_score"] >= 400]
)

string_high_confidence.to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/ppi/string_high_confidence.csv", index=False
)
