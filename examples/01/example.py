from rdkit import Chem
from neoralab_sascorer import sa_score_mol

if __name__ == "__main__":
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    print(sa_score_mol(mol))