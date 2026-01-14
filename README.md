# neoralab-sascorer

[![License](https://img.shields.io/github/license/neoralab/neoralab-sascorer)](LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/neoralab-sascorer)](https://pypi.org/project/neoralab-sascorer/)
[![PyPI](https://img.shields.io/pypi/v/neoralab-sascorer)](https://pypi.org/project/neoralab-sascorer/)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](docs/index.md)

![Synthetic accessibility scoring banner](assets/hero.png)

A small, pip-installable wrapper around RDKit's `Contrib/SA_Score` implementation for
computing the synthetic accessibility (SA) score.

📚 **Documentation:** See the MkDocs site content in [`docs/index.md`](docs/index.md).

## Features

- Simple Python API for SMILES or RDKit molecule inputs.
- CLI for quick SA score checks from the shell.
- Lightweight packaging with dynamic versioning from git tags.

## Requirements

- Python 3.8+
- RDKit (install via conda or ensure a compatible wheel is available for your platform)

## Installation

> **Dependency**: This package requires RDKit. `pip install neoralab-sascorer` expects
> an RDKit wheel to be available for your platform. If you prefer conda, you can install
> RDKit via conda and then install this package.

Editable install from this repo:

```bash
uv pip install -e .
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

## Development

### Versioning

This project uses `versioningit` to derive versions from git tags. For example, a
commit after the `0.1.0` tag will produce a version like `0.1.0.dev3`.

### Pre-commit

Install and run pre-commit hooks locally:

```bash
uv pip install pre-commit
pre-commit install
pre-commit run --all-files
```

## Attribution

This package bundles the SA_Score implementation and data originally distributed in
RDKit's `Contrib/SA_Score`, based on the method described by Ertl and Schuffenhauer
("Estimation of synthetic accessibility score of drug-like molecules based on molecular
complexity and fragment contributions", *J. Cheminformatics* 1:8, 2009).
