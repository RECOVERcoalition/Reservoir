from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def fingerprint_similarity_from_fingerprint(fp1, fp2):
    """ Compute a Tanimoto coefficient
    """

    return DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.TanimotoSimilarity)

def fingerprint_similarity_from_smiles(smiles1, smiles2):
    """ Compute the Tanimoto coefficient of the Morgan Fingerprint
    """
    try:
        fp1 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles1),2, nBits=1024)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles2),2, nBits=1024)

        return fingerprint_similarity_from_fingerprint(fp1, fp2)
    except Exception as e:
        print(f"Error {str(e)}")

if __name__ == "__main__":
    drug1 = "CCn1cc(C(=O)[O-])c(=O)c2ccc(C)nc21"
    drug2 = "CC1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23"

    print(fingerprint_similarity_from_smiles(drug1, drug2))
