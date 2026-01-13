# neoralab-sascorer

A small, pip-installable wrapper around RDKit's `Contrib/SA_Score` implementation for
computing the synthetic accessibility (SA) score.

## Installation

> **Dependency**: This package requires RDKit. `pip install neoralab-sascorer` expects
> an RDKit wheel to be available for your platform. If you prefer conda, you can install
> RDKit via conda and then install this package:
>
> ```bash
> conda install -c conda-forge rdkit
> pip install -e .
> ```

Editable install from this repo:

```bash
pip install -e .
```

The `fpscores.pkl.gz` fragment score file must be present in
`neoralab-sascorer/src/neoralab_sascorer/` (it is distributed with RDKit's
`Contrib/SA_Score` and should be copied in before building a wheel).

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

## Attribution

This package bundles the SA_Score implementation and data originally distributed in
RDKit's `Contrib/SA_Score`, based on the method described by Ertl and Schuffenhauer
("Estimation of synthetic accessibility score of drug-like molecules based on molecular
complexity and fragment contributions", *J. Cheminformatics* 1:8, 2009).
