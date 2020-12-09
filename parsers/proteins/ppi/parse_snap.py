import pandas as pd
import reservoir as rsv
from reservoir.parsers.proteins import protein_mapper
import ppi_helper

# load data
snap = pd.read_csv(rsv.RESERVOIR + "/raw/ppi/snap-msi.tsv", sep="\t")

# Node names are hgnc symbols. Update to latest version
snap["gene_1_hgnc_id"] = snap["node_1_name"].apply(protein_mapper.hgnc_normalize)
snap["gene_2_hgnc_id"] = snap["node_2_name"].apply(protein_mapper.hgnc_normalize)

# some nodes are missing a name but have the ensembl id
snap.loc[pd.isnull(snap["gene_1_hgnc_id"]), "gene_1_hgnc_id"] = snap.loc[
    pd.isnull(snap["gene_1_hgnc_id"]), "node_1"
].apply(lambda eid: protein_mapper.ensembl_to_hgnc(f"ENSG{eid:011}"))
snap.loc[pd.isnull(snap["gene_2_hgnc_id"]), "gene_2_hgnc_id"] = snap.loc[
    pd.isnull(snap["gene_2_hgnc_id"]), "node_2"
].apply(lambda eid: protein_mapper.ensembl_to_hgnc(f"ENSG{eid:011}"))

# remove duplicates
snap = ppi_helper.remove_duplicate_edges(snap[["gene_1_hgnc_id", "gene_2_hgnc_id"]])

# export
snap.to_csv(rsv.RESERVOIR + "/parsed/ppi/snap.csv", index=False)
