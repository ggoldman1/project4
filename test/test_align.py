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
            assert nw._gapB_matrix[r,c] == max(
                nw.gap_open + nw.gap_extend + nw._align_matrix[r-1, c],
                nw.gap_open + nw.gap_extend + nw._gapA_matrix[r-1, c],
                nw.gap_extend + nw._gapB_matrix[r-1, c]
            )

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

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    assert nw.align(seq3, seq4) == (17.0, 'MAVHQLIRRP', 'M---QLIRHP')

    for r in range(1, nw._align_matrix.shape[0]):
        for c in range(1, nw._align_matrix.shape[1]):
            assert nw._back[r,c] == np.argmax(
                [nw._align_matrix[r-1, c-1],
                nw._gapA_matrix[r-1, c-1],
                nw._gapB_matrix[r-1, c-1]]
            )
            assert nw._back_A[r,c] == np.argmax(
                [nw.gap_open + nw.gap_extend + nw._align_matrix[r, c-1],
                nw.gap_extend + nw._gapA_matrix[r, c-1],
                nw.gap_open + nw.gap_extend + nw._gapB_matrix[r, c-1]]
            )
            assert nw._back_B[r,c] == np.argmax(
                [nw.gap_open + nw.gap_extend + nw._align_matrix[r-1, c],
                nw.gap_open + nw.gap_extend + nw._gapA_matrix[r-1, c],
                nw.gap_extend + nw._gapB_matrix[r-1, c]]
            )




