from pymol import cmd
from collections import Counter
from pandas import DataFrame
from ..interface.pymol import bond_colors_from_array, bond_connections_from_array, select_clusters, color_selections, show_cluster
from ..interface.wordom import read_pathway_edge_frequencies
from .matrix import matrix_to_colorarray

'''
 Internal procedures
 Copyright (C) 2018  Robert Pilstål

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
'''


def draw_ciacg(cigraph, residuemap, pdb, cutoffs):
    """draw Correlated Interaction Allosteric Communication Graph (ciACG) in PyMOL

    :param cigraph: ciACG, a symmetric numpy array
    :param residuemap: OrderedDict of residue numbers to residue name mappings
    :param pdb: a pdb filename, str
    :param cutoffs: list of floats with allosteric connection strength cutoffs
    :return: list of lists with residue nodes on different cutoff levels
    """
    cmd.load(pdb)
    cmd.hide("everything")
    cmd.show("ribbon")

    # Create bindings and selections, and color them
    levels = []
    for cutoff in cutoffs:
        residues = bond_connections_from_array(cigraph, residuemap, cutoff=cutoff)
        levels.append(residues)
    selections = select_clusters(levels)
    colors = color_selections(selections)

    # Show clusters
    show_cluster(levels)
    
    return levels


def highlight_pathways(pathways, residuemap, cutoff = 0.0):
    rgb_matrix = matrix_to_colorarray(pathways)
    colored, colors = bond_colors_from_array(rgb_matrix, residuemap, cutoff = cutoff)
    return rgb_matrix, colored, colors

def process_framefiles(framefiles, residuemap):
    """Procedure to read and normalize edge counts in multiple .frames

    :param framefiles: list of strings with filenames to WORDOM .frame
                       files
    :param residuemap: dict with residue names to integer mappings
    :return: Pandas dataframe of normalized edge counts, 
             Counter of unique files processed,
             Counter of frames discovered and processed,
             Counter of unique start and endpoints discovered & proc.
    """
    files_processed = Counter()
    frames_processed = Counter()
    pathways_processed = Counter()
    frequencies = DataFrame()

    for frame in framefiles:
        files_processed[frame] += 1
        with open(frame, 'r') as infile:
            new_frequencies, new_frames, new_pathways = read_pathway_edge_frequencies(infile, residuemap)
            frequencies = frequencies.add(new_frequencies, fill_value = 0.0)
            frames_processed += new_frames
            pathways_processed += new_pathways

    frequencies = frequencies.fillna(value = 0.0)

    print("Resulting in")

    print(frequencies)

    # Normalize counts w.r.t. frames analyzed and pathways found
    # NOTE; currently counting all processed frames, while only
    #       counting those endpoints for which any shortest path
    #       were found - not considering those for which none were
    unique_frames = len(frames_processed)
    unique_pathways = len(pathways_processed)
    frequencies = frequencies.divide(unique_frames * unique_pathways)

    return frequencies, files_processed, frames_processed, pathways_processed


