from __future__ import annotations

import sys

from neoralab_sascorer.api import sa_score


def main() -> int:
    args = sys.argv[1:]
    if not args:
        return 1

    for smiles in args:
        try:
            score = sa_score(smiles)
        except ValueError as exc:
            print(f"{smiles}\tERROR: {exc}", file=sys.stderr)
            continue
        print(f"{smiles}\t{score:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
