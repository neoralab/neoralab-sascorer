from neoralab_sascorer import sa_score

if __name__ == "__main__":
    score = sa_score("CC(=O)Oc1ccccc1C(=O)O")
    print(score)