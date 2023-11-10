# -*- coding: utf-8 -*-
"""
Created on Nov 1 14:04:39 2023

@author: jaidev
"""


import streamlit as st
import numpy as np

st.title("Pairwise Sequence Alignment")
st.write("This web application performs pairwise sequence alignment using Needleman Wunsch Algorithm")

def needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty):
    # Initialize the scoring matrix
    score_matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]

    # Initialize the traceback matrix
    traceback_matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]

    # Fill the first row and column with gap penalties
    for i in range(1, len(seq1) + 1):
        score_matrix[i][0] = score_matrix[i-1][0] + gap_penalty
        traceback_matrix[i][0] = 'U'
    for j in range(1, len(seq2) + 1):
        score_matrix[0][j] = score_matrix[0][j-1] + gap_penalty
        traceback_matrix[0][j] = 'L'

    # Fill the scoring and traceback matrices
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
            if score_matrix[i][j] == match:
                traceback_matrix[i][j] = 'D'
            elif score_matrix[i][j] == delete:
                traceback_matrix[i][j] = 'U'
            else:
                traceback_matrix[i][j] = 'L'

    # Traceback to find the alignment
    alignment1 = ''
    alignment2 = ''
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 'D':
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'U':
            alignment1 = seq1[i-1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        else:
            alignment1 = '-' + alignment1
            alignment2 = seq2[j-1] + alignment2
            j -= 1

    return score_matrix[-1][-1], alignment1, alignment2


# Example usage
seq1 = st.text_input("Enter the sequence 1")
seq2 = st.text_input("Enter the sequence 2")
match_score = st.number_input("Enter the Match Score")
mismatch_score = st.number_input("Enter the Mismatch score")
gap_penalty = st.number_input("Enter the Gap Penalty")


if st.button('Submit'):
    score, alignment1, alignment2 = needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty)
    st.write("Alignment Score:", score)
    st.write("Alignment 1:", alignment1)
    st.write("Alignment 2:", alignment2)
    
    


        
    
    





