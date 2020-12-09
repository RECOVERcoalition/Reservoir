import pandas as pd
import reservoir as rsv
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover

# load the references for mapping
drug_names = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/drug_names.csv"
).rename(columns={"name": "drug_name"})
recover_to_drugbank = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_to_drugbank.csv"
)
recover_to_pubchem = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_to_pubchem.csv",
    dtype={"pubchem_cid": "str"},
)
recover_to_chembl = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_to_chembl.csv"
)

recover_drugs = pd.read_csv(
    rsv.RESEVOIR_DATA_FOLDER + "/parsed/drugs/recover_drugs.csv"
)


REFERENCES = {
    "pubchem_cid": recover_to_pubchem,
    "drugbank_id": recover_to_drugbank,
    "chembl_id": recover_to_chembl,
    "molregno": recover_to_chembl,
    "drug_name": drug_names,
    "input_smiles": recover_drugs[["smiles", "recover_id"]],
}


def map_to_recover_id_using_reference(drugs, reference):
    """ Map the drugs to a specified reference and split into mapped and unmapped
    """
    # map to reference
    drugs = drugs.merge(reference, how="left")

    # extract mapped and unmapped
    mapped_drugs = drugs.loc[~pd.isnull(drugs.recover_id)]
    unmapped_drugs = drugs.loc[pd.isnull(drugs.recover_id)]

    # remove recover_id fron unmapped ones
    del unmapped_drugs["recover_id"]

    return mapped_drugs, unmapped_drugs


def parse_smiles(smiles):
    """ Sanity check and normalization for drugs

    If a string is not provided (e.g. None, nan) None is returned back
    If a string is provided and can't be parsed it is returned back

    Args:
        smiles: a string with the smiles

    """
    if type(smiles) != type(""):
        return None

    try:
        # Copnvert to Rdkit format
        smiles = smiles.split()[0]
        mol = Chem.MolFromSmiles(smiles)

        # remove salts
        remover = SaltRemover()
        mol = remover.StripMol(mol)

        # get new smiles
        parsed_smiles = Chem.MolToSmiles(mol)
        return parsed_smiles

    except Exception as e:
        print(f"Error parsing smiles: {smiles}")

    return smiles


def map_drugs(drugs, commit_new_drugs=False):
    """ Maps a list of drugs to a Relation id

    Args:
        drugs (Pandas DataFrame): dataframe with the following possible columns
                    drug_name
                    pubchem_cid
                    drugbank_id
                    chembl_id
                    input_smiles

    Returns:
        Pandas DataFrame: The same dataframe with the Relation drug ID
    """
    # check if we have all the columns if commiting
    if commit_new_drugs and (
        "drug_name" not in drugs.columns or "input_smiles" not in drugs.columns
    ):
        raise Exception(
            "In ordered to commit new drugs you need to supply a data frame with the drug_name and input_smiles columns"
        )

    unmapped_drugs = drugs
    mapped_drugs_list = []

    # normalize the drug name for matching. the reference will be lower case
    if "drug_name" in drugs.columns:
        drugs["drug_name"] = drugs["drug_name"].str.lower()

    if "input_smiles" in drugs.columns:
        drugs["smiles"] = drugs["input_smiles"].apply(parse_smiles)

    # go through each reference. map and continue with unmapped
    for column in REFERENCES:
        if column in drugs.columns:
            # map and continue with unmapped
            mapped_drugs, unmapped_drugs = map_to_recover_id_using_reference(
                unmapped_drugs, REFERENCES[column]
            )
            mapped_drugs_list.append(mapped_drugs)

    if commit_new_drugs:
        print(unmapped_drugs)

    # return mapped and unmapped with a None in recover_id
    return pd.concat(mapped_drugs_list + [unmapped_drugs])


if __name__ == "__main__":
    pass
