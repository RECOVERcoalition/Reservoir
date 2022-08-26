import recover_data_lake as rdl
import os
import pandas as pd
import pickle
import networkx as nx
import itertools


# """
# Parse Tron
# """
# tron_folder = hgnc_data = rdl.RELATION_DATA_FOLDER + "/raw/cell_lines/tron_crawl"
# parsed_tron_files = []
# for file in os.listdir(tron_folder):
#     if not ".csv" in file:
#         continue
#
#     cell_line = file.replace(".csv", "")
#
#     try:
#         # recognize gene and keep not null ones
#         tron_file = pd.read_csv(tron_folder + "/" + file)
#         tron_file["gene_hgnc_id"] = tron_file["gene"].apply(rdl.recognize_gene_id)
#         tron_file = tron_file.loc[
#             ~tron_file["gene_hgnc_id"].isna(), ["gene_hgnc_id", "expression"]
#         ]
#
#         # some genes appear twice
#         tron_file = tron_file.sort_values(
#             "expression", ascending=False
#         ).drop_duplicates("gene_hgnc_id")
#         tron_file["cell_line"] = cell_line
#
#         parsed_tron_files.append(tron_file)
#     except:
#         print(f"Failed cell line {cell_line}")
#
#
# # combine all of the cell lines, and make the gene the column names
# parsed_tron_files = pd.concat(parsed_tron_files)
# parsed_tron_files = parsed_tron_files.set_index(["cell_line", "gene_hgnc_id"]).unstack(
#     -1
# )["expression"]
#
# # convert TPM
# parsed_tron_files = parsed_tron_files.apply(
#     lambda row: row / row.sum() * 10 ** 6, axis=1
# )
#
# parsed_tron_files = parsed_tron_files.reset_index()
# parsed_tron_files["cell_line_id"] = parsed_tron_files["cell_line"].apply(
#     map_cell_line_name
# )

from create_cell_line_db import clean_cell_line_name

with open(rdl.RECOVER_DATA_FOLDER + "/parsed/cell_lines/cd.pickle", "rb") as f:
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
    rdl.RECOVER_DATA_FOLDER + "/raw/cell_lines/E-MTAB-2770-query-results.tpms.tsv",
    sep="\t",
)
ccce["gene_hgnc_id"] = ccce["Gene Name"].apply(rdl.hgnc_normalize)
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
    rdl.RECOVER_DATA_FOLDER + "/parsed/cell_lines/cancer_cell_encyclopedia.csv",
    index=False,
)
