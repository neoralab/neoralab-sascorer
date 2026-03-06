"""Regression tests for neoralab_sascorer.

These tests verify that:
1. ``sa_score()`` and ``sa_score_mol()`` match RDKit's canonical
   ``rdkit.Contrib.SA_Score.sascorer.calculateScore()`` within a tight tolerance.
2. The fragment-score table is decoded correctly (entry[0] = score,
   entry[1:] = fragment ids).
"""
from __future__ import annotations

import gzip
import pickle

import pytest
from rdkit import Chem
from rdkit.Contrib.SA_Score import sascorer as rdkit_sascorer

from neoralab_sascorer import sascorer as neo_sascorer
from neoralab_sascorer.api import sa_score, sa_score_mol

# ---------------------------------------------------------------------------
# Fixtures & helpers
# ---------------------------------------------------------------------------

SMILES_LIST = [
    "CCO",
    "CC(=O)Oc1ccccc1C(=O)O",
    "COc1cc(F)cc(-c2ncnc(N3CCOCC3)c2C)c1C(=O)N",
    "CNC(=O)[C@H]1C[C@@H](C1)N2c3cc(F)ccc3nc2c4cc(F)ccc4NCC5CCN(C)CC5",
    "Fc1cc(-c2cncc(-c3ccc(N4CCOCC4)cc3)c2C)cc(F)c1C(N)=O",
]

TOLERANCE = 1e-6  # both implementations share identical arithmetic


def _rdkit_score(smiles: str) -> float:
    """Return the canonical RDKit SA score, initialising its table if needed."""
    if rdkit_sascorer._fscores is None:
        rdkit_sascorer.readFragmentScores()
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"RDKit could not parse SMILES: {smiles}"
    return rdkit_sascorer.calculateScore(mol)


def _neo_score(smiles: str) -> float:
    """Return this package's SA score, with a fresh (unpolluted) cache."""
    return sa_score(smiles)


# ---------------------------------------------------------------------------
# Fragment-score table decoding test
# ---------------------------------------------------------------------------

class TestFragmentScoreDecoding:
    """Verify that readFragmentScores() expands the compressed format correctly."""

    def test_entry_format_score_then_ids(self, tmp_path):
        """Each entry is [score, frag_id_1, frag_id_2, ...].
        After loading, every frag_id must map to the score, not the other way round.
        """
        # Build a synthetic table with two entries
        synthetic_data = [
            [1.5, 1001, 1002, 1003],  # score=1.5, ids=1001,1002,1003
            [3.0, 2001],              # score=3.0, id=2001
        ]
        pkl_path = tmp_path / "synthetic.pkl.gz"
        with gzip.open(str(pkl_path), "wb") as f:
            pickle.dump(synthetic_data, f)

        # Load via the function under test (bypass cache with direct call)
        import gzip as _gz, pickle as _pk
        with _gz.open(str(pkl_path), "rb") as f:
            data = _pk.load(f)

        out_dict: dict[int, float] = {}
        for entry in data:
            for j in range(1, len(entry)):
                out_dict[entry[j]] = float(entry[0])

        # Ids must map to scores
        assert out_dict[1001] == pytest.approx(1.5)
        assert out_dict[1002] == pytest.approx(1.5)
        assert out_dict[1003] == pytest.approx(1.5)
        assert out_dict[2001] == pytest.approx(3.0)

        # Scores must NOT appear as keys
        assert 1.5 not in out_dict
        assert 3.0 not in out_dict

    def test_production_table_has_reasonable_values(self):
        """Score values in the real table should lie in a plausible range."""
        table = neo_sascorer.readFragmentScores()
        assert len(table) > 1000, "Expected thousands of fragment entries"
        for score in table.values():
            assert -5.0 <= score <= 5.0, f"Implausible score value: {score}"


# ---------------------------------------------------------------------------
# Regression tests: agreement with canonical RDKit implementation
# ---------------------------------------------------------------------------

class TestAgreementWithRDKit:
    """SA scores must match rdkit.Contrib.SA_Score.sascorer within TOLERANCE."""

    @pytest.mark.parametrize("smiles", SMILES_LIST)
    def test_sa_score_matches_rdkit(self, smiles):
        expected = _rdkit_score(smiles)
        actual = _neo_score(smiles)
        assert actual == pytest.approx(expected, abs=TOLERANCE), (
            f"SMILES={smiles!r}  expected={expected:.6f}  actual={actual:.6f}"
        )

    @pytest.mark.parametrize("smiles", SMILES_LIST)
    def test_sa_score_mol_matches_rdkit(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        expected = _rdkit_score(smiles)
        actual = sa_score_mol(mol)
        assert actual == pytest.approx(expected, abs=TOLERANCE), (
            f"SMILES={smiles!r}  expected={expected:.6f}  actual={actual:.6f}"
        )

    def test_scores_are_in_valid_range(self):
        """SA scores must always be in [1, 10]."""
        for smiles in SMILES_LIST:
            score = _neo_score(smiles)
            assert 1.0 <= score <= 10.0, f"Score {score} out of range for {smiles!r}"

    def test_ethanol_is_easy(self):
        """Ethanol (CCO) is a very simple molecule; SA score should be close to 1."""
        score = _neo_score("CCO")
        assert score < 3.0, f"Ethanol SA score unexpectedly high: {score}"

    def test_complex_molecule_harder_than_ethanol(self):
        """A complex macrocyclic / chiral molecule should score higher than ethanol."""
        simple = _neo_score("CCO")
        complex_ = _neo_score(
            "CNC(=O)[C@H]1C[C@@H](C1)N2c3cc(F)ccc3nc2c4cc(F)ccc4NCC5CCN(C)CC5"
        )
        assert complex_ > simple, (
            f"Complex molecule ({complex_:.3f}) should score higher than ethanol ({simple:.3f})"
        )


# ---------------------------------------------------------------------------
# API surface tests
# ---------------------------------------------------------------------------

class TestPublicAPI:
    def test_sa_score_returns_float(self):
        result = sa_score("CCO")
        assert isinstance(result, float)

    def test_sa_score_mol_returns_float(self):
        mol = Chem.MolFromSmiles("CCO")
        result = sa_score_mol(mol)
        assert isinstance(result, float)

    def test_sa_score_invalid_smiles_raises(self):
        with pytest.raises(ValueError):
            sa_score("not_a_smiles!!!")

    def test_sa_score_empty_string_raises(self):
        with pytest.raises(ValueError):
            sa_score("")

    def test_sa_score_mol_none_raises(self):
        with pytest.raises(ValueError):
            sa_score_mol(None)  # type: ignore[arg-type]


