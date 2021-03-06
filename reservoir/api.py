import reservoir as rsv
import pandas as pd
import os
import numpy as np


def get_ppis(ppi_file, gene_hgnc_ids=set([])):
    """ Returns a specific ppi file
    """

    if type(gene_hgnc_ids) not in [type([]), type(set([]))]:
        raise Exception("Gene HGCN ids need to be provided in a list or set")

    # get ppis
    ppis = pd.read_csv(f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/ppi/{ppi_file}")

    if len(gene_hgnc_ids) > 0:
        # filtering by gene hgcn ids
        gene_hgnc_ids = set(gene_hgnc_ids)
        ppis = ppis.loc[
            ppis.apply(
                lambda row: row["gene_1_hgnc_id"] in gene_hgnc_ids
                or row["gene_2_hgnc_id"] in gene_hgnc_ids,
                axis=1,
            )
        ]

    return ppis


def available_ppis():
    """ Returns all of the available PPI files
    """

    csv_files = [
        file
        for file in os.listdir(f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/ppi/")
        if ".csv" in file
    ]

    return csv_files


def get_dtis(dti_file, recover_ids=set([])):
    """ Returns a specific dti file
    """

    if type(recover_ids) not in [type([]), type(set([]))]:
        raise Exception("Recover IDs need to be provided in a list or set")

    # get dtis
    dtis = pd.read_csv(f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/dti/{dti_file}")

    if len(recover_ids) > 0:
        # filtering by recover ids
        recover_ids = set(recover_ids)
        dtis = dtis.loc[dtis.recover_id.apply(lambda re_id: re_id in recover_ids)]

    return dtis


def available_dtis():
    """Returns all of the available DTI files
    """

    csv_files = [
        file
        for file in os.listdir(f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/dti/")
        if ".csv" in file
    ]

    return csv_files


from reservoir.parsers.drugs import molecule_mapper


def get_drugs(filter=set([])):
    """Returns the drugs in the Recover database with names, smiles and relation ids
    """

    if type(filter) not in [type([]), type(set([]))]:
        raise Exception("Recover IDs need to be provided in a list or set")

    recover_drugs = pd.read_csv(
        f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/drugs/recover_drugs.csv"
    )

    if len(filter) > 0:
        # filtering by recover ids
        filter = set(filter)
        recover_drugs = recover_drugs.loc[
            recover_drugs.recover_id.apply(lambda re_id: re_id in filter)
        ]

    return recover_drugs


def get_proteins():
    """Returns all of the available human proteins in the Relation database
    """

    return pd.read_csv(
        f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/proteins/human_proteins.csv"
    )


def get_covid_proteins():
    """ Returns all of the covid proteins
    """

    return pd.read_csv(
        f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/proteins/covid_proteins.csv"
    )


def get_gene_ontology_evidence_categories():
    """ Get the types of evidence for the GO annotations
    """

    go = pd.read_csv(f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/proteins/gene_ontology.csv")
    return list(go["evidence_category"].unique())


def get_gene_ontology(
    gene_hgnc_ids=set([]),
    go_ids=set([]),
    evidence_category=None,
    aspect=None,
    max_distance_from_root=None,
):
    """
    Returns gene ontology annotations. These can be filtered by genes or by aspect

    Args:
        gene_hgnc_ids: a list or set of gene ids
        go_ids: a list of go ids
        aspect: "cellular_component" or "biological_process" or "molecular_function"
        max_distance_from_root: int indicate max hops from root
        evidence_category: string filter for the type of evidence for the annotation. sugest using "experimental evidence"
    """

    if type(gene_hgnc_ids) not in [type([]), type(set([]))]:
        raise Exception("Gene HGCN ids need to be provided in a list or set")

    if type(go_ids) not in [type([]), type(set([]))]:
        raise Exception("Gene HGCN ids need to be provided in a list or set")

    go = pd.read_csv(f"{rsv.RESERVOIR_DATA_FOLDER}/parsed/proteins/gene_ontology.csv")

    # get rid of roots and also non-connected components (-1)
    go = go.loc[go.distance_from_root > 0]

    # apply gene filter if present
    if len(gene_hgnc_ids) > 0:
        go = go.loc[
            go.gene_hgnc_id.apply(lambda gene_hgnc_id: gene_hgnc_id in gene_hgnc_ids)
        ]

    # apply go term filter if present
    if len(go_ids) > 0:
        go = go.loc[go.go_id.isin(go_ids)]

    # apply aspect filter if present
    if aspect != None:
        go = go.loc[go.aspect == aspect]

    # apply max_distance_from_root id present
    if not max_distance_from_root is None:
        go = go.loc[go.distance_from_root <= max_distance_from_root]

    # apply evidence category filter if present
    if not evidence_category is None:
        go = go.loc[go.evidence_category == evidence_category]

    return go


def get_gene_onotology_embeddings(distance_from_root=2, evidence_category=None):
    """
    Create embeddings for proteins based on wether a protein has a label or not

    Args:
        distance_from_root: int indicating at which distance from the root we should make the embeddings
        evidence_category: string filter for the type of evidence for the annotation. sugest using "experimental evidence"
    """

    # get human proteins
    human_proteins = rsv.get_proteins()

    # extract the required go annotations
    go = rsv.get_gene_ontology(evidence_category=evidence_category)
    go = go.loc[go.distance_from_root == distance_from_root]

    # creates embedding by checking each unique term if it's present
    def create_go_embedding(terms_for_gene, unique_terms):
        terms_for_gene = set(terms_for_gene["go_id"])
        return np.array([1 if term in terms_for_gene else 0 for term in unique_terms])

    # identify all terms and then create embedding from annotations
    unique_terms = go["go_id"].unique()
    embeddings = (
        go.groupby("gene_hgnc_id")
        .apply(lambda rows: create_go_embedding(rows, unique_terms))
        .reset_index()
        .rename(columns={0: "embedding"})
    )

    # some proteins do not have annotations. create embeddings with all 0s
    missing_proteins = set(human_proteins["gene_hgnc_id"]).difference(
        set(embeddings["gene_hgnc_id"])
    )
    missing_proteins_embeddings = pd.DataFrame(
        [
            {"gene_hgnc_id": gene, "embedding": np.array([0] * len(unique_terms))}
            for gene in missing_proteins
        ]
    )

    return pd.concat([embeddings, missing_proteins_embeddings])


def get_specific_drug_combo_blocks(
    mono_row_measurements=None,
    mono_col_measurements=None,
    combo_measurements=None,
    study_name=None,
    cell_line_name=None,
):
    """ Used to pre-select combo blocks by metadata.
    """
    summary_data = pd.read_csv(
        rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/summary_data.csv",
        low_memory=False,
    )

    if not mono_row_measurements is None:
        summary_data = summary_data.loc[
            summary_data["mono_row_measurements"] == mono_row_measurements
        ]

    if not mono_col_measurements is None:
        summary_data = summary_data.loc[
            summary_data["mono_col_measurements"] == mono_col_measurements
        ]

    if not combo_measurements is None:
        summary_data = summary_data.loc[
            summary_data["combo_measurements"] == combo_measurements
        ]

    if not study_name is None:
        summary_data = summary_data.loc[summary_data["study_name"] == study_name]

    if not cell_line_name is None:
        summary_data = summary_data.loc[
            summary_data["cell_line_name"] == cell_line_name
        ]

    return summary_data


def get_drug_combo_data_monos(block_ids=None):
    mono_data = pd.read_json(
        rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/mono_therapy.json"
    )

    if not block_ids is None:
        mono_data = mono_data.loc[mono_data.block_id.isin(set(block_ids))]

    return mono_data


def get_drug_combo_data_combos(block_ids=None):
    combo_data = pd.read_json(
        rsv.RESERVOIR_DATA_FOLDER + "/parsed/drug_combos/combos.json"
    )

    if not block_ids is None:
        combo_data = combo_data.loc[combo_data.block_id.isin(set(block_ids))]

    return combo_data
