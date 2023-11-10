# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 14:04:39 2023

@author: jaide
"""


import streamlit as st
import numpy as np

st.title("Pairwise Sequence Alignment")
st.write("This web application performs pairwise sequence alignment using Smith Waterman Algorithm. This a local sequence alignment algorithm.")

def smith_waterman(seq1, seq2, match_score, mismatch_score, gap_penalty):
    # Initialize the scoring matrix
    score_matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    
    # Initialize the traceback matrix
    traceback_matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    
    # Initialize the maximum score and its position
    max_score = 0
    max_pos = (0, 0)
    
    # Fill the scoring and traceback matrices
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(0, match, delete, insert)
            
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)
                
                if score_matrix[i][j] == match:
                    traceback_matrix[i][j] = 1
                elif score_matrix[i][j] == delete:
                    traceback_matrix[i][j] = 2
                else:
                    traceback_matrix[i][j] = 3
    
    # Traceback to find the alignment
    alignment_seq1 = ""
    alignment_seq2 = ""
    i, j = max_pos
    
    while score_matrix[i][j] > 0:
        if traceback_matrix[i][j] == 1:
            alignment_seq1 = seq1[i-1] + alignment_seq1
            alignment_seq2 = seq2[j-1] + alignment_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 2:
            alignment_seq1 = seq1[i-1] + alignment_seq1
            alignment_seq2 = "-" + alignment_seq2
            i -= 1
        else:
            alignment_seq1 = "-" + alignment_seq1
            alignment_seq2 = seq2[j-1] + alignment_seq2
            j -= 1
    
    return max_score, alignment_seq1, alignment_seq2


# Example usage
seq1 = st.text_input("Enter the sequence 1")
seq2 = st.text_input("Enter the sequence 2")
match_score = st.number_input("Enter the Match Score")
mismatch_score = st.number_input("Enter the Mismatch score")
gap_penalty = st.number_input("Enter the Gap Penalty")


if st.button('Submit'):
    score, alignment1, alignment2 = smith_waterman(seq1, seq2, match_score, mismatch_score, gap_penalty)
    st.write("Alignment Score:", score)
    st.write("Alignment 1:", alignment1)
    st.write("Alignment 2:", alignment2)
    
    


        
    
    





