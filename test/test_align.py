# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
nw.sub_dict = {k:1 if k[0] == k[1] else -1 for k in nw.sub_dict.keys()}
nw.gap_open = -3
    nw.align(seq1, seq2)
    assert nw.align(seq1, seq2) == (4.0, 'MYQR', 'M-QR'), "alignment score and/or alignment is wrong"

    for r in range(1, nw._align_matrix.shape[0]):
        for c in range(1, nw._align_matrix.shape[1]):
            assert nw._align_matrix[r,c] == max(
                nw._align_matrix[r-1, c-1],
                nw._gapA_matrix[r-1, c-1],
                nw._gapB_matrix[r-1, c-1]
            ) + nw.sub_dict[(seq1[r-1], seq2[c-1])], "alignment matrix has incorrect values"
            assert nw._gapA_matrix[r,c] == max(
                nw.gap_open + nw.gap_extend + nw._align_matrix[r, c-1],
                nw.gap_extend + nw._gapA_matrix[r, c-1],
                nw.gap_open + nw.gap_extend + nw._gapB_matrix[r, c-1]
            )




    assert np.all(align_mat)

    gapA = nw._gapA_matrix == [[-11., -np.inf, -np.inf, -np.inf],
                               [-12., -13., -6., -7.],
                               [-13., -14., -15., -7.],
                               [-14., -15., -16., -12.],
                               [-15., -16., -17., -17.]]
    assert np.all(gapA)

    gapB = nw._gapB_matrix == [[-11., -12., -13., -14.],
                               [-np.inf, -13., -14., -15.],
                               [-np.inf,  -6., -15., -16.],
                               [-np.inf,  -7.,  -7., -17.],
                               [-np.inf,  -8.,  -8.,  -6.]]
    assert np.all(gapB)


def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    pass




