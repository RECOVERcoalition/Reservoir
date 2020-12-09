import pandas as pd
import reservoir_data_lake as rsv
from reservoir_data_lake.parsers.proteins import protein_mapper
import ppi_helper

# read in data and map protein id
huri_original = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/ppi/HuRI.tsv",
    sep="\t",
    names=["protein_1_ensembl", "protein_2_ensembl"],
)

# normalize gene ids and keep only valid ones
huri_original["gene_1_hgnc_id"] = huri_original["protein_1_ensembl"].apply(
    protein_mapper.ensembl_to_hgnc
)
huri_original["gene_2_hgnc_id"] = huri_original["protein_2_ensembl"].apply(
    protein_mapper.ensembl_to_hgnc
)
huri_original = huri_original.loc[
    ~pd.isnull(huri_original.gene_1_hgnc_id) & ~pd.isnull(huri_original.gene_2_hgnc_id)
]

# remove duplicate edges
huri_original = ppi_helper.remove_duplicate_edges(huri_original)

huri_original.to_csv(rsv.RESERVOIR_DATA_FOLDER + "/parsed/ppi/huri.csv", index=False)
