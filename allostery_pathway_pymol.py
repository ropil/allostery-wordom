#!/usr/bin/env python3
if __name__ == "__main__" and __package__ is None:
    __package__ = "allostery-wordom"

import pymol
from pymol import cmd

from .internal.procedure import draw_ciacg, highlight_pathways, process_framefiles
from .internal.matrix import align_dataframes, matrix_from_pandas_dataframe

import numpy
import matplotlib.pyplot as plt
from pandas import DataFrame

import pickle
'''
 Display PSNPath on a ciACG in an interactive PyMOL session
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


# Version and license information
def get_version_str():
    return "\n".join([
        "allostery_pathway_pymol  Copyright (C) 2018  Robert Pilstål;",
        "This program comes with ABSOLUTELY NO WARRANTY.",
        "This is free software, and you are welcome to redistribute it",
        "under certain conditions; See the supplied Apache License,",
        "Version 2.0 for the specific language governing permissions",
        "and limitations under the License.",
        "    http://www.apache.org/licenses/LICENSE-2.0"
    ])


# Library functions
#def get_cutoffs(cutoffs):
    #return [float(cutoff) for cutoff in cutoffs.split(',')]


# Main; for callable scripts
def main():
    from argparse import ArgumentParser
    from sys import argv
    parser = ArgumentParser(
        description= "Open and display a WORDOM PSN and cross correlation " +
                     "analysis as an correlated interaction Allosteric " +
                     "Communication Graph (ciACG) in PyMOL."
    )
    parser.add_argument(
        "-c",
        nargs='*',
        default=[0.0],
        metavar="float",
        help="Cutoff to use for showing connections in graph, provide" +
        " white space separated list for a series of cutoffs."
        ", default=0.0")
    parser.add_argument(
        "-pdb", nargs=1, metavar="PDBfile", help="PDB file to draw")
    parser.add_argument("-plot", action="store_true", default=False, help="Plot ciACG value distribution")
    parser.add_argument(
        "-acg", nargs=1, metavar="ACGfile", help="ACG file to read (.frm)")
    parser.add_argument(
        "-rmp", nargs=1, metavar="RMPfile", help="ResidueMap file to read (.rmp)")
    parser.add_argument(
        "-pml", nargs='*', metavar="PMLfile", default=None, help="PyMOL scripts to run with cmd.run(), before coloring of bonds.")
    parser.add_argument(
        "-frames", nargs='*', metavar="FRAMEfile", help="WORDOM .frame files to process")
    arguments = parser.parse_args(argv[1:])

    # Finish pymol launch
    pymol.finish_launching(['pymol'])

    # Set variables here
    pdb = arguments.pdb[0]
    cutoffs = [float(c) for c in arguments.c]
    acg = arguments.acg[0]
    rmp = arguments.rmp[0]
    pml = arguments.pml
    frames = arguments.frames

    with open(acg, 'rb') as infile:
        cigraph_table = pickle.load(infile)

    with open(rmp, 'rb') as infile:
        residuemap = pickle.load(infile)

    frequencies, files_processed, frames_processed, pathways_processed = process_framefiles(frames, residuemap)

    print(frequencies)

    print("{} pathways found in {} frames from {} files".format(len(pathways_processed), len(frames_processed), len(files_processed)))

    # Align tables
    frequencies_aligned, cigraph_table_aligned = align_dataframes(frequencies, cigraph_table, fill_value = 0.0)

    pathways = matrix_from_pandas_dataframe(frequencies_aligned)

    print("Aligned frequencies")
    print(frequencies_aligned)
    print("Aligned cigraph")
    print(cigraph_table_aligned)
    print("Residuemap")
    print(residuemap)

    # Draw the loaded ciACG
    cigraph = matrix_from_pandas_dataframe(cigraph_table_aligned)
    levels = draw_ciacg(cigraph, residuemap, pdb, cutoffs)

    # Run scripts prior to coloring of bonds
    if pml is not None:
        for script in pml:
            cmd.run(script)

    # Highlight the pathways
    rgb_matrix, colored, colors = highlight_pathways(pathways, residuemap)

if __name__ == '__main__':
    main()
