import reservoir as rsv
import pandas as pd
import pickle


def create_alias_dictionary(hgnc_data):
    """ Create a dictionary mapping each historical gene id to the current gene ids
    """
    alias_dictionary = {}
    for i, row in hgnc_data.iterrows():
        if type(row["Alias symbol"]) == type(""):
            alias_dictionary[row["Alias symbol"]] = row["Approved symbol"]
        if type(row["Previous symbol"]) == type(""):
            alias_dictionary[row["Previous symbol"]] = row["Approved symbol"]

    return alias_dictionary


hgnc_data = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/proteins/hgnc_output.tab", sep="\t"
)

# select protein-producing genes that have a uniprot id
hgnc_data = hgnc_data.loc[
    (hgnc_data.Status == "Approved")
    & (hgnc_data["Locus type"] == "gene with protein product")
    & ~pd.isnull(hgnc_data["UniProt accession"])
]

# make all symbols uppercase to avoid issues
hgnc_data["Approved symbol"] = hgnc_data["Approved symbol"].apply(
    lambda symbol: symbol.upper() if type(symbol) == type("") else None
)
hgnc_data["Alias symbol"] = hgnc_data["Alias symbol"].apply(
    lambda symbol: symbol.upper() if type(symbol) == type("") else None
)
hgnc_data["Previous symbol"] = hgnc_data["Previous symbol"].apply(
    lambda symbol: symbol.upper() if type(symbol) == type("") else None
)

# create dictionary for mapping previous symbols to the current one
alias_dictionary = create_alias_dictionary(hgnc_data)
with open(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/gene_alias_dictionary.pickle", "wb"
) as f:
    pickle.dump(alias_dictionary, f)

# add protein sequence
hgnc_data = hgnc_data[
    [
        "HGNC ID",
        "Status",
        "Approved symbol",
        "Approved name",
        "UniProt accession",
        "Ensembl gene ID",
    ]
].drop_duplicates()
uniprot_data = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/raw/proteins/uniprot_output.tab", sep="\t"
)
hgnc_data = hgnc_data.merge(
    uniprot_data[["Entry", "Sequence"]], left_on="UniProt accession", right_on="Entry"
)

# export genes
hgnc_data = hgnc_data.rename(
    columns={
        "Approved symbol": "gene_hgnc_id",
        "Approved name": "gene_name",
        "UniProt accession": "uniprot_id",
        "Ensembl gene ID": "gene_ensembl_id",
        "Sequence": "protein_sequence",
    }
)
hgnc_data[
    ["gene_hgnc_id", "gene_name", "gene_ensembl_id", "uniprot_id", "protein_sequence"]
].to_csv(rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/human_proteins.csv", index=False)
