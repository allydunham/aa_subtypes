#!/usr/bin/env python3
"""
Module containing functions to calculate structural environment profiles of AAs
"""
# TODO Tests to make sure these are producing correct profiles

from collections import defaultdict
import numpy as np
from Bio.SeqUtils import seq1
from Bio.Alphabet.IUPAC import protein as protein_alphabet

AA_INDEX = {aa: i for i, aa in enumerate(protein_alphabet.letters)}
NUM_AAS = len(protein_alphabet.letters)

def k_nearest_residues(residues, k=10, distance_matrix=None):
    """
    Yields chemical environments parameterised by the make up of the k nearest AAs.
    Hetero atoms are included so must be dropped separately if desired.

    residues: list of residues to consider
    k: count the k nearest residues
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order. Calculated if not supplied

    yields: chemical environment profiles (np.array)
    """
    if k >= len(residues):
        raise ValueError('k >= number of residues')

    if distance_matrix is None:
        distance_matrix = residue_distance_matrix(residues)

    for res_index in range(len(residues)):
        dists = distance_matrix[res_index,]

        non_self = np.ones_like(dists, dtype=bool)
        non_self[res_index] = False

        nearest_k = [residues[i] for i in np.argpartition(dists[non_self], k)[:k]]

        counts = defaultdict(lambda: 0)
        for i in nearest_k:
            counts[seq1(i.get_resname())] += 1

        yield np.array([counts[aa] for aa in protein_alphabet.letters])

# 10A selected as default because used in Bagley & Altman 1995
def within_distance(residues, max_dist=10, distance_matrix=None):
    """
    Yeilds chemical environments parameterised as the residues within max_dist
    angstroms. Hetero atoms are included so must be dropped separately if desired.

    residues: list of residues to consider
    max_dist: maximum distance to count within (in Angstroms)
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order. Calculated if not supplied

    yields: chemical environment profiles (np.array)
    """

    if distance_matrix is None:
        distance_matrix = residue_distance_matrix(residues)

    for res_index in range(len(residues)):
        dists = distance_matrix[res_index,]

        res_within_dist = [residues[i] for i in np.argwhere(dists < max_dist)[:, 0]
                           if not i == res_index]

        counts = defaultdict(lambda: 0)
        for i in res_within_dist:
            counts[seq1(i.get_resname())] += 1

        yield np.array([counts[aa] for aa in protein_alphabet.letters])

def distance_to_nearest(residues, distance_matrix=None):
    """
    Yeilds chemical environments parameterised as the distance to the nearest residue
    of each type. Hetero atoms are included so must be dropped separately if desired.

    residues: list of residues to consider
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order. Calculated if not supplied

    yields: chemical environment profiles (np.array)
    """

    if distance_matrix is None:
        distance_matrix = residue_distance_matrix(residues)

    residue_indices = [np.array([seq1(r.get_resname()) == aa for r in residues]) for
                       aa in protein_alphabet.letters]

    for res_index in range(len(residues)):
        dists = distance_matrix[res_index,]

        non_self = np.ones_like(dists, dtype=bool)
        non_self[res_index] = False

        yield np.array([min(dists[aa & non_self]) if any(aa & non_self) else np.nan for
                        aa in residue_indices])

def residue_distance_matrix(residues, ref_atom='CA'):
    """
    Generate a distance matrix from an iterable of Bio.PDB residues.
    There is no checking whether the specified atom makes sense (i.e. using C-beta
    with Gly will fail)

    residues: iterable of Bio.PDB residues to get distances from
    ref_atom: atom to measure distances from

    returns: distance matrix, rows and columns in order of residues (np.array)
    """

    dist = np.zeros((len(residues), len(residues)))

    for i, res1 in enumerate(residues):
        for j, res2 in enumerate(residues[i + 1:]):
            dist[i, j + i + 1] = dist[j + i + 1, i] = res1[ref_atom] - res2[ref_atom]

    return dist
