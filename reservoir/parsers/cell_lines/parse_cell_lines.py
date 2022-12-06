import reservoir as rsv
import os
import pandas as pd
import pickle
import networkx as nx
import itertools


from create_cell_line_db import clean_cell_line_name

with open(rsv.RESERVOIR_DATA_FOLDER + "/parsed/cell_lines/cd.pickle", "rb") as f:
    cd = pickle.load(f)


def map_cell_line_name(cell_line_name):
    clean_name = clean_cell_line_name(cell_line_name)
    if clean_name in cd:
        return cd[clean_name]

    return None


"""
EBI RNA-seq of cancer cell lines
"""
# load and recognize gene ids. keep ones with valid ids
ccce = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/cell_lines/E-MTAB-2770-query-results.tpms.tsv",
    sep="\t",
)
ccce["gene_hgnc_id"] = ccce["Gene Name"].apply(rsv.hgnc_normalize)
ccce = ccce.loc[~ccce["gene_hgnc_id"].isna()]
ccce = ccce[ccce.columns[2:]]

# put the cell line as the index
ccce = ccce.set_index("gene_hgnc_id").T.reset_index()
ccce["index"] = ccce["index"].apply(lambda index: index.split(", ")[1])
ccce = ccce.rename(columns={"index": "cell_line"})

# some values are missing, fill them in
ccce = ccce.fillna(0)

# normalize to TPM
ccce = (
    ccce.set_index("cell_line")
    .apply(lambda row: row / row.sum() * 10 ** 6, axis=1)
    .reset_index()
)

# set internal id
ccce["cell_line_id"] = ccce["cell_line"].apply(map_cell_line_name)
ccce = ccce.loc[~ccce["cell_line_id"].isna()]

ccce.to_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/cell_lines/cancer_cell_encyclopedia.csv",
    index=False,
)
