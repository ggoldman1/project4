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

    assert np.all(nw._align_matrix == [[  0., -np.inf, -np.inf, -np.inf],
                                       [-np.inf,   5., -11., -13.],
                                       [-np.inf, -12.,   4.,  -8.],
                                       [-np.inf, -12.,  -1.,   5.],
                                       [-np.inf, -14.,  -6.,   4.]])

    assert np.all(nw._gapA_matrix == [[-10., -np.inf, -np.inf, -np.inf],
                                      [-11., -12.,  -6.,  -7.],
                                      [-12., -13., -14.,  -7.],
                                      [-13., -14., -15., -12.],
                                      [-14., -15., -16., -17.]])

    assert np.all(nw._gapB_matrix == [[-10., -11., -12., -13.],
                                      [-np.inf, -12., -13., -14.],
                                      [-np.inf,  -6., -14., -15.],
                                      [-np.inf,  -7.,  -7., -16.],
                                      [-np.inf,  -8.,  -8.,  -6.]])

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

    assert np.all(nw._back == [[-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf],
                               [-np.inf,   0.,   2.,   2.,   2.,   2.,   2.,   2.],
                               [-np.inf,   1.,   0.,   1.,   1.,   1.,   1.,   1.],
                               [-np.inf,   1.,   2.,   0.,   0.,   0.,   0.,   1.],
                               [-np.inf,   1.,   2.,   2.,   0.,   0.,   1.,   1.],
                               [-np.inf,   1.,   2.,   0.,   2.,   0.,   0.,   0.],
                               [-np.inf,   1.,   2.,   0.,   2.,   0.,   0.,   0.],
                               [-np.inf,   1.,   2.,   2.,   0.,   0.,   2.,   0.],
                               [-np.inf,   1.,   2.,   2.,   0.,   0.,   1.,   1.],
                               [-np.inf,   1.,   2.,   0.,   2.,   2.,   0.,   1.],
                               [-np.inf,   1.,   2.,   0.,   2.,   2.,   0.,   0.]])

    assert np.all(nw._back_A == [[-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf],
                                 [-np.inf,   1.,   0.,   1.,   1.,   1.,   1.,   1.],
                                 [-np.inf,   1.,   1.,   0.,   1.,   1.,   1.,   1.],
                                 [-np.inf,   1.,   1.,   1.,   0.,   1.,   1.,   1.],
                                 [-np.inf,   1.,   1.,   1.,   2.,   0.,   1.,   0.],
                                 [-np.inf,   1.,   1.,   0.,   1.,   1.,   0.,   1.],
                                 [-np.inf,   1.,   1.,   1.,   0.,   1.,   1.,   0.],
                                 [-np.inf,   1.,   1.,   1.,   0.,   0.,   1.,   1.],
                                 [-np.inf,   1.,   1.,   0.,   2.,   2.,   0.,   1.],
                                 [-np.inf,   1.,   1.,   0.,   2.,   2.,   0.,   0.],
                                 [-np.inf,   1.,   1.,   1.,   2.,   2.,   2.,   2.]])

    assert np.all(nw._back_B == [[-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf],
                                 [-np.inf,   2.,   2.,   2.,   2.,   2.,   2.,   2.],
                                 [-np.inf,   0.,   2.,   2.,   2.,   2.,   2.,   2.],
                                 [-np.inf,   2.,   0.,   2.,   2.,   2.,   2.,   2.],
                                 [-np.inf,   2.,   2.,   0.,   0.,   1.,   1.,   1.],
                                 [-np.inf,   2.,   2.,   2.,   0.,   0.,   0.,   0.],
                                 [-np.inf,   2.,   2.,   2.,   2.,   0.,   2.,   0.],
                                 [-np.inf,   2.,   2.,   2.,   2.,   2.,   0.,   2.],
                                 [-np.inf,   2.,   2.,   2.,   0.,   2.,   2.,   2.],
                                 [-np.inf,   2.,   2.,   2.,   2.,   0.,   1.,   1.],
                                 [-np.inf,   2.,   2.,   2.,   2.,   2.,   0.,   1.]])




