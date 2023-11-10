# -*- coding: utf-8 -*-
"""
Created on Nov 10 2023

@author: jaidev
"""
import streamlit as st
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import time


def calculate_distance(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    return len(best_alignment[0]) - best_alignment[2]

def upgma(sequences):
    n = len(sequences)
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            distances[i][j] = calculate_distance(sequences[i], sequences[j])
    linkage_matrix = linkage(distances, method='average') 
    fig = plt.figure()
    plt.xlabel('Sequences')
    plt.ylabel('Distance')
    plt.title('UPGMA Tree')
    dendrogram(linkage_matrix, labels=organism_names, leaf_rotation=90)
    st.pyplot(fig)
    
st.title("Phylogenetics")
st.write("This web application is for the phylogenetic analysis of multiple sequences. It forms the distance tree between various sequences belonging to different organisms based on the UPGMA Algorithm (Unweighted Pair Group Method with Arithmetic Mean). To use this application, input the sequences and the labels for the sequences of different organisms in Space Separated format")
sequ = st.text_input("Enter the sequences in space separated format")
labe = st.text_input("Enter the labels witth space between the names of different organisms")
sequences = sequ.split(" ")
organism_names = labe.split(" ")
if st.button("Submit"):
    with st.spinner('Running the UPGMA to form the distance tree...'):
        time.sleep(5)
    st.success('Done!')
    upgma(sequences)
