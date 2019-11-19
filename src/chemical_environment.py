#!/usr/bin/env python3
"""
Module containing functions to calculate structural environment profiles of AAs
"""

# TODO Tests to make sure these are producing correct profiles

import sys
from collections import defaultdict
import numpy as np
from Bio.SeqUtils import seq1#, seq3
from Bio.Alphabet.IUPAC import protein as protein_alphabet

AA_INDEX = {aa: i for i, aa in enumerate(protein_alphabet.letters)}
NUM_AAS = len(protein_alphabet.letters)

class ChemicalEnvironment:
    """
    Object determining chemical environments for residues in a Bio.PDB model
    """
    def __init__(self, model, pdb_id=None, ref_atom='CA'):
        self.pdb_id = model.get_full_id()[0] if pdb_id is None else pdb_id
        self.ref_atom = ref_atom
        self.model = model

        self.residues = None
        self.residue_groups = None
        self.distance_matrices = None
        self._dist_matrix_generated = False

        self.k_nearest_profiles = {}
        self.within_distance_profiles = {}

    def select_residues(self, groups='all'):
        """
        Select groups of residues to profile within. Accepted values for groups:
            - 'all': group all residues from all chains [default]
            - 'chains': separate each chain into a group
            - list of lists of chain letters to group
        More complex groupings should be assigned manually to self.residues with
        self.residue_groups assigned as a list of appropriate names

        groups: groups of residues to process together. Structured as a list of dicts
                with each dict having a key for each included chain, mapping to an iteralble of
                residue indecies to include or None for the whole chain.
                If groups = None all residues will be grouped and the special value 'chains'
                will separate each chain
        drop_hetero: logical determining if hetero atoms will be filtered
        """
        self.residue_groups = groups
        if groups is None:
            self.residues = [[r for c in self.model for r in c.get_residues()]]

        elif groups == 'chains':
            self.residues = [list(c.get_residues()) for c in self.model]

        else:
            self.residues = [[r for chain in group for r in self.model[chain].get_residues()] for
                             group in groups]

        self.residues = [drop_hetero_atoms(r) for r in self.residues]

        self._dist_matrix_generated = False

    def _make_distance_matrices(self):
        """
        Generate distance matrices for each group of residues
        """
        self.distance_matrices = [residue_distance_matrix(r, ref_atom=self.ref_atom) for
                                  r in self.residues]
        self._dist_matrix_generated = True

    # TODO check there are more than k residues in set
    def calc_k_nearest_profiles(self, k):
        """
        Calculate k-nearest AAs chemical environment profile
        """
        if self.residues is None:
            raise ValueError('No residues selected')

        if not self._dist_matrix_generated or self.distance_matrices is None:
            self._make_distance_matrices()

        self.k_nearest_profiles[k] = [[k_nearest_residues_profile(i, res, mat, k) for
                                       i in range(len(res))] for
                                      res, mat in zip(self.residues, self.distance_matrices)]

    def calc_within_distance_profiles(self, max_dist):
        """
        Calculate 'AAs within max_dist angstroms' chemical environment profile
        """
        if self.residues is None:
            raise ValueError('No residues selected')

        if not self._dist_matrix_generated or self.distance_matrices is None:
            self._make_distance_matrices()

        profs = [[within_distance_profile(i, res, mat, max_dist) for
                  i in range(len(res))] for
                 res, mat in zip(self.residues, self.distance_matrices)]
        self.within_distance_profiles[max_dist] = profs

    def aa_profiles(self, offset=0):
        """
        Generator yielding each AA chemical environment

        offset: int to correct PDB positions by, either as a dict with values per
                chain or an integer
        """
        if isinstance(offset, int):
            _offset = defaultdict(lambda: offset)
        else:
            _offset = offset

        chain_ids = [x.get_id() for x in self.model.get_list()]
        for i, residues in enumerate(self.residues):
            # generate group name
            if self.residue_groups == 'all':
                group = 'all'
            elif self.residue_groups == 'chains':
                group = f'chain_{chain_ids[i]}'
            elif isinstance(self.residue_groups[i], str):
                group = self.residue_groups[i]
            else:
                group = (f"chain{'s' if len(self.residue_groups[i]) > 1 else ''}"
                         f"_{'_'.join(self.residue_groups[i])}")

            for j, res in enumerate(residues):
                res_id = res.get_full_id()
                yield {'pdb_id': self.pdb_id, 'chain': res_id[2], 'aa': seq1(res.get_resname()),
                       'position': int(res_id[3][1]) + _offset[res_id[2]], 'group': group,
                       **{f'nearest_{k}': v[i][j] for k, v in self.k_nearest_profiles.items()},
                       **{f'within_{k}': v[i][j] for k, v in self.within_distance_profiles.items()}}

    def write_profile_table(self, file_obj, offset=0):
        """
        Write all profiles in tsv format

        file_obj: file like object to write to
        offset: int to correct PDB positions by, either as a dict with values per
                chain or an integer. Passed to self.aa_profiles()
        """
        col_names = ['pdb_id', 'chain', 'group', 'position', 'aa']
        prof_col_start = len(col_names)
        col_names.extend([f'nearest_{k}' for k in self.k_nearest_profiles])
        col_names.extend([f'within_{k}' for k in self.within_distance_profiles])

        print(*col_names, sep='\t', file=file_obj)
        for aa_prof in self.aa_profiles(offset=offset):
            for key in col_names[prof_col_start:]:
                aa_prof[key] = ','.join(map(str, aa_prof[key]))

            print(*[aa_prof[k] for k in col_names], sep='\t', file=file_obj)

def k_nearest_residues_profile(res_index, residues, distance_matrix=None, k=10):
    """
    Calculate chemical environment of a residue as the make up of the k nearest AAs.
    Hetero atoms are included so must be dropped separately if desired.

    res_index: residue to process
    residues: list of all residues to consider
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order
    k: count the k nearest residues

    returns: structural profile of the residue (np.array)
    """
    if k >= len(residues):
        print('Warning: k >= number of residues. Returning NA', sys.stderr)
        return np.repeat(np.nan, NUM_AAS)

    if distance_matrix is None:
        distance_matrix = residue_distance_matrix(residues)

    dists = distance_matrix[res_index,]
    dists[res_index] = dists.max() + 1 # don't want the residue itself

    nearest_k = [residues[i] for i in np.argpartition(dists, k)[:k]]

    counts = defaultdict(lambda: 0)
    for i in nearest_k:
        counts[seq1(i.get_resname())] += 1

    return np.array([counts[aa] for aa in protein_alphabet.letters])

# 10A selected as default because used in Bagley & Altman 1995
def within_distance_profile(res_index, residues, distance_matrix=None, max_dist=10):
    """
    Calculate chemical environment of a residue as the residues within max_dist
    angstroms of it. Hetero atoms are included so must be dropped separately if desired.

    res_index: residue to process
    residues: list of all residues to consider
    distance_matrix: numpy matrix of distances between residues, with rows/columns in
                     that order
    max_dist: maximum distance to count within

    returns: structural profile of the residue (np.array)
    """

    if distance_matrix is None:
        distance_matrix = residue_distance_matrix(residues)

    dists = distance_matrix[res_index,]
    dists[res_index] = max_dist + 1 # ignore residue itself

    res_within_dist = [residues[i] for i in np.argwhere(dists < max_dist)[:, 0]]

    counts = defaultdict(lambda: 0)
    for i in res_within_dist:
        counts[seq1(i.get_resname())] += 1

    return np.array([counts[aa] for aa in protein_alphabet.letters])

def get_profiles(residues, profile_func=k_nearest_residues_profile):
    """
    Calculate chemical environment profiles for a list of residues

    residues: list of Bio.PDB residues to retrieve profiles for
    profile_func: function to calculate profiles with sig f(res_ind, res_list, dist_matrix)
    drop_hetero: filter out hetero atoms (bool)

    returns: profiles (list of np.array)
    """
    residues = drop_hetero_atoms(residues)

    distance_matrix = residue_distance_matrix(residues)

    return [profile_func(i, residues, distance_matrix) for i in range(len(residues))]

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

def drop_hetero_atoms(residues):
    """
    Drop hetero atoms.

    residues: iterable of Bio.PDB residues

    returns: residues without hetero atoms (list)
    """
    return [r for r in residues if r.id[0] == ' ']
