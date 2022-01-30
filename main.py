# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    score_against_hs = {}
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    for seq, species in zip([gg_seq, mm_seq, br_seq, tt_seq], ["Gallus gallus", "Mus musculus",
                                                               "Balaencipes rex", "Tursiops truncatus"]):
        score_against_hs[species] = nw.align(hs_seq, seq)[0]
    # sort species so highest scorer goes first 
    score_against_hs_sorted = {k: v for k, v in sorted(score_against_hs.items(),
                                                       key=lambda item: item[1], reverse=True)}
    print("Species in order of similarity to human:")
    for spec in score_against_hs_sorted.keys():
        print(f"\t{spec}")

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print(f"\nAlignment scores against human:")
    for spec in score_against_hs_sorted.keys():
        print(f"\t{spec}: {score_against_hs_sorted[spec]}")

if __name__ == "__main__":
    main()
