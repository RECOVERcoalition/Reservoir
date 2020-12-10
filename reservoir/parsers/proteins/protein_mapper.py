import pandas as pd
import pickle
import reservoir as rsv

# load human proteins and alias mapper
human_proteins = pd.read_csv(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/human_proteins.csv"
)
approved_human_proteins = set(human_proteins["gene_hgnc_id"])
with open(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/gene_alias_dictionary.pickle", "rb"
) as f:
    gene_alias_dictionary = pickle.load(f)

# create mappings from different type of IDs to HGNC
uniprot_to_hgnc_map = human_proteins.set_index("uniprot_id").to_dict()["gene_hgnc_id"]
ensembl_to_hgnc_map = human_proteins.set_index("gene_ensembl_id").to_dict()[
    "gene_hgnc_id"
]
with open(
    rsv.RESERVOIR_DATA_FOLDER + "/parsed/proteins/covid_protein_mapper.pickle", "rb"
) as f:
    covid_protein_mapper = pickle.load(f)


def map_covid_protein(protein_id):
    """ Maps one of the multiple covid protein identifiers to the one we use to index
    """

    protein_id = protein_id.lower()

    if protein_id in covid_protein_mapper:
        return covid_protein_mapper[protein_id]

    return None


def hgnc_normalize(hgnc_symbol):
    """
    Maps a previous hgnc symbol to the current. If it doesn't match to anything return None

    Args:
        hgnc_symbol
    """
    hgnc_symbol = str(hgnc_symbol).upper()

    if hgnc_symbol in gene_alias_dictionary:
        hgnc_symbol = gene_alias_dictionary[hgnc_symbol]

    if hgnc_symbol in approved_human_proteins:
        return hgnc_symbol
    else:
        return None


def uniprot_to_hgnc(uniprot_id):
    """ Maps a uniprot id to the HGCN gene symbol
    """

    if uniprot_id in uniprot_to_hgnc_map:
        return uniprot_to_hgnc_map[uniprot_id]
    else:
        return None


def ensembl_to_hgnc(ensembl_id):
    """ Maps an ensembl gene id to the HGCN gene symbol
    """

    if ensembl_id in ensembl_to_hgnc_map:
        return ensembl_to_hgnc_map[ensembl_id]
    else:
        return None
