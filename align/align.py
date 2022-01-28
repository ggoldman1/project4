# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        # TODO: Fill in the Needleman-Wunsch Algorithm below
        to perform global sequence alignment of seqA and seqB
        and return a tuple with the following format
        (alignment score, seqA alignment, seqB alignment)
        Also, write up a docstring for this function using the
        _read_sub_matrix as an example.
        Don't forget to comment your code!
        """
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        # number each `back` matrix for ease of use
        self._back_dict = {0: self._back, 1: self._back_A, 2: self._back_B}

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        #### NEEDLEMAN-WUNSCH IMPLEMENTATION ####
        # Begin with base cases
        self._align_matrix[0,0] = 0
        self._gapA_matrix[:,0] = [self.gap_open + (i+1)*self.gap_extend for i in range(len(seqA) + 1)]
        self._gapB_matrix[0,:] = [self.gap_open + (i+1)*self.gap_extend for i in range(len(seqB) + 1)]

        # Now, fill in the "interior" of each matrix
        for row in range(1, self._align_matrix.shape[0]):
            for col in range(1, self._align_matrix.shape[1]):

                max_align = self._max_step_align(row, col)
                self._align_matrix[row, col] = self.sub_dict[(seqA[row-1], seqB[col-1])] + max_align[0]
                self._back = max_align[1], max_align[2] # store the matrix and the element

                max_gapA = self._max_step_gapA(row, col)
                self._gapA_matrix[row, col] = max_gapA[0]
                self._back_A = max_gapA[1], max_gapA[2]

                max_gapB = self._max_step_gapB(row, col)
                self._gapB_matrix[row, col] = max_gapB[0]
                self._back_B = max_gapB[1], max_gapB[2]

        return self._backtrace()

    def _get_max(self, align: int, gapA: int, gapB, i: int, j: int) -> Tuple:
        """
        Helper method for filling in matrices. Gets max of three elements, and returns the matrix (and coordinates)
        containing the max value.

        Parameters:
            align: int
                Alignment option
            gapA: int
                GapA option
            gapB: int
                GapB option
            i: int
                Current row
            j: int
                Current column

        Return:
            Tuple
                Maximum element from above
                Which element contained max element {0: align, 1: gapA, 2: gapB}
                Which entry in the element contained the max element

        """
        if align > gapA and align > gapB:
            return (align, 0, (i-1, j-1))
        elif gapA > align and gapA > gapB:
            return (gapA, 1, (i-1, j-1))
        return (gapB, 2, (i-1, j-1))

    def _max_step_align(self, i: int, j: int) -> Tuple:
        """
        Helper method for filling in M[i,j], ie calculates
            -----------
            | M[i-1, j-1]
        max | gapA[i-1, j-1]
            | gapB[i-1, j-1]
            -----------
        It also returns which matrix and which element within the matrix contained the maximum element, for back tracing.

        Parameters:
            i: int
                Current row of iteration
            j: int
                Current column of iteration

        Returns:
            Tuple
                Maximum element from above
                Which element contained max element {0: align, 1: gapA, 2: gapB}
                Which entry in the element contained the max element
        """
        align = self._align_matrix[i-1, j-1]
        gapA = self._gapA_matrix[i-1, j-1]
        gapB = self._gapB_matrix[i-1, j-1]

        return self._get_max(align, gapA, gapB, i, j)

    def _max_step_gapA(self, i: int, j: int) -> Tuple:
        """
        Helper method for filling in gapA[i,j], ie calculates
            -----------
            | gap_start + gap_extent + M[i, j-1]
        max | gap_extend + gapA[i, j-1]
            | gap_start + gap_extend + gapB[i, j-1]
            -----------
        It also returns which matrix and which element within the matrix contained the maximum element, for back tracing.

        Parameters:
            i: int
                Current row of iteration
            j: int
                Current column of iteration

        Returns:
            Tuple
                Maximum element from above
                Which element contained max element {0: align, 1: gapA, 2: gapB}
                Which entry in the element contained the max element
        """
        align = self.gap_open + self.gap_extend + self._align_matrix[i, j-1]
        gapA = self.gap_extend + self._gapA_matrix[i, j-1]
        gapB = self.gap_open + self.gap_extend + self._gapB_matrix[i, j-1]

        return self._get_max(align, gapA, gapB, i, j)

    def _max_step_gapB(self, i: int, j: int) -> Tuple:
        """
        Helper method for filling in gapB[i,j], ie calculates
            -----------
            | gap_start + gap_extent + M[i-1, j]
        max | gap_start + gap_extend + gapA[i-1, j]
            | gap_extend + gapB[i-1, j]
            -----------
        It also returns which matrix and which element within the matrix contained the maximum element, for back tracing.

        Parameters:
            i: int
                Current row of iteration
            j: int
                Current column of iteration

        Returns:
            Tuple
                Maximum element from above
                Which element contained max element {0: align, 1: gapA, 2: gapB}
                Which entry in the element contained the max element
        """
        align = self.gap_open + self.gap_extend + self._align_matrix[i-1, j]
        gapA = self.gap_open + self.gap_extend + self._gapA_matrix[i-1, j]
        gapB = self.gap_extend + self._gapB_matrix[i-1, j]

        return self._get_max(align, gapA, gapB, i, j)


    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """
        # Implement this method based upon the heuristic chosen in the align method above.
        pass


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
