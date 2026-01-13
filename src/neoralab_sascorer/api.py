from __future__ import annotations

from rdkit import Chem

from neoralab_sascorer import sascorer


def sa_score(smiles: str) -> float:
    """Return the synthetic accessibility score for a SMILES string."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("SMILES must be a non-empty string")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    return sascorer.calculateScore(mol)


def sa_score_mol(mol: Chem.Mol) -> float:
    """Return the synthetic accessibility score for an RDKit Mol."""
    if mol is None:
        raise ValueError("Mol must not be None")

    return sascorer.calculateScore(mol)
