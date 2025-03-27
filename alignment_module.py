import numpy as np

class SequenceAligner:
    def __init__(self, match_score=2, mismatch_score=-1, gap_penalty=-1):
        """
        Initialize alignment parameters
        
        Args:
            match_score (int): Score for matching nucleotides
            mismatch_score (int): Penalty for mismatched nucleotides
            gap_penalty (int): Penalty for inserting gaps
        """
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty

    def needleman_wunsch(self, seq1, seq2):
        """
        Global sequence alignment using Needleman-Wunsch algorithm
        
        Args:
            seq1 (str): First DNA/protein sequence
            seq2 (str): Second DNA/protein sequence
        
        Returns:
            dict: Alignment results with score, aligned sequences, etc.
        """
        # Create scoring matrix
        m, n = len(seq1), len(seq2)
        score_matrix = np.zeros((m+1, n+1), dtype=int)
        
        # Initialize first row and column
        for i in range(m+1):
            score_matrix[i][0] = i * self.gap_penalty
        for j in range(n+1):
            score_matrix[0][j] = j * self.gap_penalty
        
        # Fill scoring matrix
        for i in range(1, m+1):
            for j in range(1, n+1):
                match = score_matrix[i-1][j-1] + (
                    self.match_score if seq1[i-1] == seq2[j-1] 
                    else self.mismatch_score
                )
                delete = score_matrix[i-1][j] + self.gap_penalty
                insert = score_matrix[i][j-1] + self.gap_penalty
                
                score_matrix[i][j] = max(match, delete, insert)
        
        # Traceback to find alignment
        align1, align2 = [], []
        i, j = m, n
        while i > 0 and j > 0:
            score_current = score_matrix[i][j]
            score_diagonal = score_matrix[i-1][j-1]
            score_up = score_matrix[i][j-1]
            score_left = score_matrix[i-1][j]
            
            if score_current == score_diagonal + (
                self.match_score if seq1[i-1] == seq2[j-1] 
                else self.mismatch_score
            ):
                align1.insert(0, seq1[i-1])
                align2.insert(0, seq2[j-1])
                i -= 1
                j -= 1
            elif score_current == score_left + self.gap_penalty:
                align1.insert(0, seq1[i-1])
                align2.insert(0, '-')
                i -= 1
            else:
                align1.insert(0, '-')
                align2.insert(0, seq2[j-1])
                j -= 1
        
        # Handle remaining sequences
        while i > 0:
            align1.insert(0, seq1[i-1])
            align2.insert(0, '-')
            i -= 1
        while j > 0:
            align1.insert(0, '-')
            align2.insert(0, seq2[j-1])
            j -= 1
        
        return {
            'score': score_matrix[m][n],
            'aligned_seq1': ''.join(align1),
            'aligned_seq2': ''.join(align2),
            'similarity_percentage': self._calculate_similarity(
                ''.join(align1), ''.join(align2)
            )
        }

    def smith_waterman(self, seq1, seq2):
        """
        Local sequence alignment using Smith-Waterman algorithm
        
        Args:
            seq1 (str): First DNA/protein sequence
            seq2 (str): Second DNA/protein sequence
        
        Returns:
            dict: Local alignment results
        """
        m, n = len(seq1), len(seq2)
        score_matrix = np.zeros((m+1, n+1), dtype=int)
        
        # Variables to track maximum score location
        max_score = 0
        max_i, max_j = 0, 0
        
        # Fill scoring matrix
        for i in range(1, m+1):
            for j in range(1, n+1):
                match = score_matrix[i-1][j-1] + (
                    self.match_score if seq1[i-1] == seq2[j-1] 
                    else self.mismatch_score
                )
                delete = score_matrix[i-1][j] + self.gap_penalty
                insert = score_matrix[i][j-1] + self.gap_penalty
                
                score_matrix[i][j] = max(0, match, delete, insert)
                
                # Track maximum score
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_i, max_j = i, j
        
        # Traceback
        align1, align2 = [], []
        i, j = max_i, max_j
        while i > 0 and j > 0 and score_matrix[i][j] > 0:
            score_current = score_matrix[i][j]
            score_diagonal = score_matrix[i-1][j-1]
            score_up = score_matrix[i][j-1]
            score_left = score_matrix[i-1][j]
            
            if score_current == score_diagonal + (
                self.match_score if seq1[i-1] == seq2[j-1] 
                else self.mismatch_score
            ):
                align1.insert(0, seq1[i-1])
                align2.insert(0, seq2[j-1])
                i -= 1
                j -= 1
            elif score_current == score_left + self.gap_penalty:
                align1.insert(0, seq1[i-1])
                align2.insert(0, '-')
                i -= 1
            else:
                align1.insert(0, '-')
                align2.insert(0, seq2[j-1])
                j -= 1
        
        return {
            'score': max_score,
            'aligned_seq1': ''.join(align1),
            'aligned_seq2': ''.join(align2),
            'similarity_percentage': self._calculate_similarity(
                ''.join(align1), ''.join(align2)
            )
        }

    def _calculate_similarity(self, seq1, seq2):
        """
        Calculate similarity percentage between aligned sequences
        
        Args:
            seq1 (str): First aligned sequence
            seq2 (str): Second aligned sequence
        
        Returns:
            float: Similarity percentage
        """
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
        total_length = len(seq1)
        return (matches / total_length) * 100 if total_length > 0 else 0