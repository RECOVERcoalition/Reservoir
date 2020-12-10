import networkx as nx
import pandas as pd


def remove_duplicate_edges(data):
    """ Removes duplicates edges in a PPI coming from A-B and B-A

    """
    # add to nx graph and export
    data = data.drop_duplicates()
    data_nx = nx.from_edgelist(
        data.apply(lambda row: (row["gene_1_hgnc_id"], row["gene_2_hgnc_id"]), axis=1)
    )
    return pd.DataFrame(
        [{"gene_1_hgnc_id": p1, "gene_2_hgnc_id": p2} for (p1, p2) in data_nx.edges]
    )
