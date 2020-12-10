import reservoir as rsv
import pandas as pd
from multiprocessing import Pool
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors


"""
1. Connect drugs based on shared parsed smiles
"""


def parse_smiles(smiles):
    """ Sanity check and normalization for drugs """
    try:
        # Remove salts
        smiles = smiles.split()[0]
        mol = Chem.MolFromSmiles(smiles)
        remover = SaltRemover()
        mol = remover.StripMol(mol)
        parsed_smiles = Chem.MolToSmiles(mol)
        return parsed_smiles

    except Exception as e:
        pass

    return smiles


def add_parsed_smiles(all_drugs_with_smiles):
    """ Add the parsed smiles in multiprocessing fashion
    """

    # identify all unique smiles and run parsed on them
    unique_smiles = set(all_drugs_with_smiles["canonical_smiles"].unique())
    with Pool(7) as p:
        parsed_smiles = p.map(parse_smiles, unique_smiles)
    canonical_to_parsed = dict(zip(unique_smiles, parsed_smiles))

    # add parsed version
    all_drugs_with_smiles["parsed_smiles"] = all_drugs_with_smiles[
        "canonical_smiles"
    ].apply(lambda cs: canonical_to_parsed[cs])

    return all_drugs_with_smiles


all_drugs = pd.read_csv(rsv.RESEVOIR_DATA_FOLDER + "/raw/dti/chembl/all_drugs.csv")

# remove drugs with no smiles and parse smiles with our internal parser
all_drugs_with_smiles = all_drugs.loc[~pd.isnull(all_drugs["canonical_smiles"])]
all_drugs_with_smiles = add_parsed_smiles(all_drugs_with_smiles)

# find drugs that share smiles strings
all_drugs_with_smiles_grouped = all_drugs_with_smiles.groupby("parsed_smiles").agg(
    {"molregno": set}
)
all_drugs_with_smiles_grouped["count"] = all_drugs_with_smiles_grouped[
    "molregno"
].apply(lambda m: len(m))
all_drugs_with_smiles_grouped = all_drugs_with_smiles_grouped.loc[
    all_drugs_with_smiles_grouped["count"] > 1
]


"""
2. Connect drugs based on salt
"""


def map_to_parent(row):
    if row["parent_molregno"] != row["molregno"]:
        return row["parent_molregno"]

    if row["active_molregno"] != row["molregno"]:
        return row["active_molregno"]

    return None


# read in drugs salts hierarchy. The parent becomes either the active or the parent
drug_salt = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/raw/dti/chembl/drug_salt.csv", error_bad_lines=False
)
drug_salt["parent"] = drug_salt.apply(map_to_parent, axis=1)


"""
3. Put all connections on a graph
"""
import networkx as nx
import itertools

# create graph and connect using chembl molregno
mol_graph = nx.Graph()
mol_graph.add_nodes_from(list(all_drugs["molregno"]))

# add edges from smiles connections
for i, row in all_drugs_with_smiles_grouped.iterrows():
    for (mol1, mol2) in itertools.combinations(row["molregno"], 2):
        mol_graph.add_edge(mol1, mol2)

# add edges from smiles connections
for i, row in drug_salt.iterrows():
    mol_graph.add_edge(row["molregno"], row["parent"])


"""
4. Identify the representative of each unique drug
"""


def compute_mw(smiles):
    try:
        smiles = smiles.split()[0]
        mol = Chem.MolFromSmiles(smiles)
        return Descriptors.ExactMolWt(mol)
    except:
        return None


# each connected component in the graph becomes one recover id
rel_id_mapper = {}
rel_counter = 1
for cc in nx.connected_components(mol_graph):
    rel_id = f"RE-MOL-{rel_counter:010}"
    for id in cc:
        rel_id_mapper[id] = rel_id
    rel_counter += 1

# add recover ID to all drugs
all_drugs_with_smiles["recover_id"] = all_drugs_with_smiles["molregno"].apply(
    lambda molregno: rel_id_mapper[molregno]
)

# add mw to molecules. in a group of molecules, the one with lowest mw becomes the representative
unique_smiles = all_drugs_with_smiles["parsed_smiles"].unique()
with Pool(7) as p:
    mws = p.map(compute_mw, unique_smiles)
smiles_to_mw = dict(zip(unique_smiles, mws))
all_drugs_with_smiles["mw"] = all_drugs_with_smiles["parsed_smiles"].apply(
    lambda smiles: smiles_to_mw[smiles]
)
all_drugs_with_smiles = all_drugs_with_smiles.sort_values("mw")


def combine_by_recover_id(molecules):
    """ Combine multiple molecules into 1
    """

    # The top molecule, the one with lowest mw, becomes the representative
    parent_molecule = molecules.iloc[0]
    if type(parent_molecule["pref_name"]) == type(""):
        name = parent_molecule["pref_name"]
    else:
        name = parent_molecule["recover_id"]

    return pd.DataFrame(
        [
            {
                "name": name,
                "smiles": parent_molecule["parsed_smiles"],
                "mw": parent_molecule["mw"],
                "max_phase": molecules["max_phase"].max(),
            }
        ]
    )


"""
5. Identify mappings to other dbs
"""

# combine drugs and export
recover_drugs = (
    all_drugs_with_smiles.groupby("recover_id")
    .apply(combine_by_recover_id)
    .reset_index()
)
recover_drugs = recover_drugs[["recover_id", "name", "smiles", "mw", "max_phase"]]
recover_drugs.loc[
    recover_drugs["name"].apply(lambda name: "RE-MOL" not in name), "name"
] = recover_drugs.loc[
    recover_drugs["name"].apply(lambda name: "RE-MOL" not in name), "name"
].apply(
    lambda name: name.lower()
)
recover_drugs.to_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_drugs.csv", index=False
)

# Chembl mapping
recover_to_chembl = all_drugs_with_smiles[
    ["recover_id", "molregno", "chembl_id"]
].drop_duplicates()
recover_to_chembl.to_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_to_chembl.csv", index=False
)

# Other mappings
chembl_to_drugbank = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/raw/dti/chembl/chembl_to_drugbank.txt", sep="\t"
).rename(columns={"From src:'1'": "chembl_id", "To src:'2'": "drugbank_id"})
chembl_to_lincs = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/raw/dti/chembl/chembl_to_lincs.txt", sep="\t"
).rename(columns={"From src:'1'": "chembl_id", "To src:'25'": "lincs_id"})
chembl_to_pubchem = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/raw/dti/chembl/chembl_to_pubchem.txt", sep="\t"
).rename(columns={"From src:'1'": "chembl_id", "To src:'22'": "pubchem_cid"})

recover_to_drugbank = recover_to_chembl.merge(chembl_to_drugbank)[
    ["recover_id", "drugbank_id"]
]
recover_to_lincs = recover_to_chembl.merge(chembl_to_lincs)[["recover_id", "lincs_id"]]
recover_to_pubchem = recover_to_chembl.merge(chembl_to_pubchem)[
    ["recover_id", "pubchem_cid"]
]


recover_to_drugbank.to_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_to_drugbank.csv", index=False
)
recover_to_lincs.to_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_to_lincs.csv", index=False
)
recover_to_pubchem.to_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_to_pubchem.csv", index=False
)

"""
6. Create name to recover id mapping
"""

# create name to recover id mapping
drug_names = recover_drugs.loc[
    recover_drugs["name"].apply(lambda name: "RE-MOL" not in name),
    ["recover_id", "name"],
].drop_duplicates()

# add synonyms
drug_synonyms = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/raw/dti/chembl/drug_synonyms.csv",
    error_bad_lines=False,
)
drug_synonyms["synonyms"] = drug_synonyms["synonyms"].apply(lambda name: name.lower())
drug_synonyms["recover_id"] = drug_synonyms["molregno"].apply(
    lambda molregno: rel_id_mapper[molregno]
)
drug_synonyms = drug_synonyms.rename(columns={"synonyms": "name"})[
    ["recover_id", "name"]
]

# combine drug names
drug_names = pd.concat([drug_names, drug_synonyms]).drop_duplicates()

# create drug names
drug_names.to_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/drug_names.csv", index=False
)

aws_api.push_folder(rsv.RESEVOIR_DATA_FOLDER + "/parsed", "parsed")
