#
# Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
"""Synthetic accessibility score implementation."""
from __future__ import annotations

import gzip
import math
import os
import pickle
from importlib import resources
from typing import Dict, Tuple

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors

_fscores: Dict[int, float] | None = None


def _resolve_scores_path(name: str) -> str:
    local_dir = os.path.dirname(__file__)
    return os.path.join(local_dir, name)


def readFragmentScores(name: str = "fpscores") -> Dict[int, float]:
    """Load fragment scores from a pickled, gzipped file.

    Each entry in the pickle data has the shape ``[score, frag_id_1, frag_id_2, ...]``.
    ``entry[0]`` is the score; every subsequent element is a fragment id that maps
    to that score.  This matches the canonical RDKit SA_Score implementation.
    """
    global _fscores
    if _fscores is not None:
        return _fscores

    filename = name if name.endswith(".pkl.gz") else f"{name}.pkl.gz"

    try:
        with resources.open_binary(__package__, filename) as handle:
            with gzip.GzipFile(fileobj=handle) as gz:
                data = pickle.load(gz)
    except FileNotFoundError:
        fname = _resolve_scores_path(filename)
        with gzip.open(fname, "rb") as f:
            data = pickle.load(f)

    # entry[0] = score value; entry[1:] = fragment ids that share that score.
    out_dict: Dict[int, float] = {}
    for entry in data:
        for j in range(1, len(entry)):
            out_dict[entry[j]] = float(entry[0])

    _fscores = out_dict
    return _fscores


def numBridgeheadsAndSpiro(mol: Chem.Mol, ri=None) -> Tuple[int, int]:
    """Return the number of bridgehead and spiro atoms."""
    if ri is None:
        ri = mol.GetRingInfo()
    n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    n_bridgeheads = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return n_bridgeheads, n_spiro


def calculateScore(mol: Chem.Mol) -> float:
    """Calculate the synthetic accessibility score for a molecule."""
    if mol is None:
        raise ValueError("Mol is required")

    fscores = readFragmentScores()
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    fp = mfpgen.GetSparseCountFingerprint(mol)
    fps = fp.GetNonzeroElements()
    score1 = 0.0
    nf = 0
    for bit_id, v in fps.items():
        nf += v
        score1 += fscores.get(bit_id, -4.0) * v
    if nf:
        score1 /= float(nf)

    n_atoms = mol.GetNumAtoms()
    n_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    n_bridgeheads, n_spiro = numBridgeheadsAndSpiro(mol)
    n_macrocycles = sum(1 for ring in mol.GetRingInfo().AtomRings() if len(ring) > 8)

    size_penalty = n_atoms**1.005 - n_atoms
    stereo_penalty = math.log10(n_chiral_centers + 1.0)
    spiro_penalty = math.log10(n_spiro + 1.0)
    bridge_penalty = math.log10(n_bridgeheads + 1.0)
    macrocycle_penalty = math.log10(2.0) if n_macrocycles > 0 else 0.0

    score2 = 0.0 - size_penalty - stereo_penalty - spiro_penalty - bridge_penalty - macrocycle_penalty

    score3 = 0.0
    if n_atoms > len(fps):
        score3 = 0.5 * math.log(float(n_atoms) / len(fps))

    sascore = score1 + score2 + score3

    min_sa = -4.0
    max_sa = 2.5
    sascore = 11.0 - (sascore - min_sa + 1.0) / (max_sa - min_sa) * 9.0

    if sascore > 8.0:
        sascore = 8.0 + math.log(sascore + 1.0 - 9.0)
    if sascore > 10.0:
        sascore = 10.0
    elif sascore < 1.0:
        sascore = 1.0

    return sascore


if __name__ == "__main__":
    import sys

    mols = [Chem.MolFromSmiles(x) for x in sys.argv[1:]]
    for m in mols:
        if m is None:
            print("Invalid SMILES")
            continue
        print(calculateScore(m))
