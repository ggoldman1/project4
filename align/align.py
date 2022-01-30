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
        Implementation of Needleman-Wunsch algorithm, which generates an alignment between two sequences using dynamic
        programming.

        Parameters:
            seqA: str
                This is the first sequence we consider when aligning
            seqB: str
                This is the second sequence we consider when aligning
        """
        # ensure both sequences only contain letters in the alphabet -- in this case, the letters in the
        # substitution matrix
        alphabet = set(np.array([[k[0], k[1]] for k in nw.sub_dict.keys()]).flatten())
        for l in seqA + seqB:
            if l not in alphabet:
                raise ValueError(f"Your sequence has a letter ({l}) that is not found in the substitution matrix")

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
        self._gapA_matrix[:,0] = [self.gap_open + i*self.gap_extend for i in range(len(seqA) + 1)]
        self._gapB_matrix[0,:] = [self.gap_open + i*self.gap_extend for i in range(len(seqB) + 1)]

        # Now, fill in the "interior" of each matrix
        for row in range(1, self._align_matrix.shape[0]):
            for col in range(1, self._align_matrix.shape[1]):

                max_align = self._max_step_align(row, col) # get the max of the three options
                self._align_matrix[row, col] = self.sub_dict[(seqA[row-1], seqB[col-1])] + max_align[0]
                self._back[row, col] = max_align[1] # point to the matrix containing max element

                max_gapA = self._max_step_gapA(row, col) # get the max of the three options
                self._gapA_matrix[row, col] = max_gapA[0]
                self._back_A[row, col] = max_gapA[1] # point to the matrix containing max element

                max_gapB = self._max_step_gapB(row, col) # get the max of the three options
                self._gapB_matrix[row, col] = max_gapB[0]
                self._back_B[row, col] = max_gapB[1] # point to the matrix containing max element

        return self._backtrace()

    def _max_step_align(self, i: int, j: int) -> Tuple:
        """
        Helper method for filling in M[i,j], ie calculates
            -----------
            | M[i-1, j-1]
        max | gapA[i-1, j-1]
            | gapB[i-1, j-1]
            -----------
        It also returns which matrix contained the maximum element, for back tracing.

        Parameters:
            i: int
                Current row of iteration
            j: int
                Current column of iteration

        Returns:
            Tuple
                Maximum element from above
                Which matrix contained max element {0: align, 1: gapA, 2: gapB}
        """
        align = self._align_matrix[i-1, j-1]
        gapA = self._gapA_matrix[i-1, j-1]
        gapB = self._gapB_matrix[i-1, j-1]
        opts = [align, gapA, gapB]

        return opts[np.argmax(opts)], np.argmax(opts)

    def _max_step_gapA(self, i: int, j: int) -> Tuple:
        """
        Helper method for filling in gapA[i,j], ie calculates
            -----------
            | gap_start + gap_extend + M[i, j-1]
        max | gap_extend + gapA[i, j-1]
            | gap_start + gap_extend + gapB[i, j-1]
            -----------
        It also returns which matrix contained the maximum element, for back tracing.
        The elements we consdier here are ``lef'', ie (i,j-1) with respect to given (i,j).


        Parameters:
            i: int
                Current row of iteration
            j: int
                Current column of iteration

        Returns:
            Tuple
                Maximum element from above
                Which matrix contained max element {0: align, 1: gapA, 2: gapB}
        """
        align = self.gap_open + self.gap_extend + self._align_matrix[i, j-1]
        gapA = self.gap_extend + self._gapA_matrix[i, j-1]
        gapB = self.gap_open + self.gap_extend + self._gapB_matrix[i, j-1]
        opts = [align, gapA, gapB]

        return opts[np.argmax(opts)], np.argmax(opts)

    def _max_step_gapB(self, i: int, j: int) -> Tuple:
        """
        Helper method for filling in gapB[i,j], ie calculates
            -----------
            | gap_start + gap_extend + M[i, j-1]
        max | gap_start + gap_extend + gapA[i-1, j-1]
            | gap_extend + gapB[i, j-1]
            -----------
        It also returns which matrix contained the maximum element, for back tracing.
        The elements we consider here are ``up'', ie (i-1, j) with respect to given (i,j)

        Parameters:
            i: int
                Current row of iteration
            j: int
                Current column of iteration

        Returns:
            Tuple
                Maximum element from above
                Which matrix contained max element {0: align, 1: gapA, 2: gapB}
        """
        align = self.gap_open + self.gap_extend + self._align_matrix[i-1, j]
        gapA = self.gap_open + self.gap_extend + self._gapA_matrix[i-1, j]
        gapB = self.gap_extend + self._gapB_matrix[i-1, j]
        opts = [align, gapA, gapB]

        return opts[np.argmax(opts)], np.argmax(opts)

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        Constructs the alignment of the two sequences by traversing through the three pointer matrices, each iteration
        deciding whehter to align the sequences or add a gap character.

        Returns:
            Tuple
                Alignment score
                Alignment representation of sequence A
                Alignment representation of sequence B
        """
        # Implement this method based upon the heuristic chosen in the align method above.
        curr_row = len(self._seqA)
        curr_col = len(self._seqB)

        # get bottom right value in each matrix
        align_val = self._align_matrix[curr_row, curr_col]
        gapA_val = self._gapA_matrix[curr_row, curr_col]
        gapB_val = self._gapB_matrix[curr_row, curr_col]

        # find out which matrix had the max value, keep track of that matrix
        # tie breaks in order of preference are (a) align, (b), gapA, then (c) gapB. 
        curr_mat = np.argmax([align_val, gapA_val, gapB_val])

        # iterate over the characters, deciding when to align or gap until finished
        while curr_row > 0 or curr_col > 0:
            # curr_row and curr_col will decrease depending on how we gap/align
            curr_mat, curr_row, curr_col = self._backtrace_sequence_appender(curr_mat, curr_row, curr_col)

        return self._align_matrix[-1, -1], self.seqA_align, self.seqB_align

    def _backtrace_sequence_appender(self, back_mat: int, i: int, j: int) -> Tuple:
        """
        Updates the sequence alignments depending on the matrix (ie append or gap). 
        
        Parameters:
            back_mat: int
                Either a 0 (align), 1 (gap A), or 2 (gap B)
                    0 : align, move up diagonally
                    1: gap in sequence A, consume a letter of B, move left (i, j-1)
                    2: gap in sequence B, consume a letter of A, move up (i-1, j)

            i: int 
                Current spot in seq A
            j: int 
                Current spot in seq B

        Returns:
            Tuple
                Next back mat (either 0, 1, or 2), and updated i and j.
        """
        if back_mat == 0: # if we are aligning at this step
            self.seqA_align = self._seqA[i-1] + self.seqA_align # add a letter
            self.seqB_align = self._seqB[j-1] + self.seqB_align # add a letter
            # return the previous matrix (stored in align_mat), move according to rule
            nextrow_nextcol = self._back_trace_rule(self._back_dict[0][i, j], i, j)
            return self._back_dict[0][i,j], nextrow_nextcol[0], nextrow_nextcol[1]

        elif back_mat == 1: # if we are gapping sequence A
            # if we are gapping at sequence A
            self.seqA_align = '-' + self.seqA_align # add a gap character to sequence A
            self.seqB_align = self._seqB[j-1] + self.seqB_align # add letter to sequence B
            # return the previous matrix pointer (stored in gapA), move according to rule
            nextrow_nextcol = self._back_trace_rule(self._back_dict[0][i, j], i, j)
            return self._back_dict[0][i, j], nextrow_nextcol[0], nextrow_nextcol[1]

        # if we are gapping at sequence B
        self.seqA_align = self._seqA[i - 1] + self.seqA_align  # add letter to sequence A
        self.seqB_align = '-' + self.seqB_align  # add a gap character to sequence be
        # return the previous matrix pointer (stored in gapB), move according to rule
        nextrow_nextcol = self._back_trace_rule(self._back_dict[0][i, j], i, j)
        return self._back_dict[0][i, j], nextrow_nextcol[0], nextrow_nextcol[1]


    def _back_trace_rule(self, back_mat: int, i: int, j: int) -> Tuple[int, int]:
        """
        Determine how to move through the alignment matrix based on the next decision (align or gap).
            0: move up diagonally
            1: move left (i, j-1)
            2: move up (i-1, j)

        Parameters:
            back_mat: int
                What is the next move (0, 1, or 2?)
            i: int
                Current row (ie letter in seqA)
            j: int
                Current col (ie letter in seqB)

        Returns:
            Tuple:
                Updated i, j based on back_mat
        """
        if back_mat == 0:
            return i-1, j-1
        elif back_mat == 1:
            return i, j-1
        else:
            return i-1, j


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
