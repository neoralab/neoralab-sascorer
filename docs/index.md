# neoralab-sascorer

Synthetic accessibility scoring from RDKit Contrib SA_Score, packaged as a small
Python library with a CLI for quick checks.

## Why use this package?

- Compute the SA score for SMILES strings or RDKit molecules.
- Reuse the original RDKit `Contrib/SA_Score` logic in a pip-installable package.
- Batch-score molecules from the command line.

## Installation

> **Dependency:** RDKit is required. If you use conda, install RDKit first and then
> install this package with pip/uv.

```bash
uv pip install neoralab-sascorer
```

To contribute locally:

```bash
uv pip install -e .
```

## Usage

### Python API

```python
from neoralab_sascorer import sa_score

score = sa_score("CC(=O)Oc1ccccc1C(=O)O")
print(score)
```

```python
from rdkit import Chem
from neoralab_sascorer import sa_score_mol

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
print(sa_score_mol(mol))
```

### CLI

```bash
neoralab-sascorer "CC(=O)Oc1ccccc1C(=O)O"
```

Output format:

```
<SMILES>\t<score>
```

## Project resources

- [Repository README](https://github.com/neoralab/neoralab-sascorer/blob/main/README.md)
