#!/usr/bin/env python3
if __name__ == "__main__" and __package__ is None:
    __package__ = "allostery-wordom"

import pymol
from pymol import cmd

from .interface.pymol import (bond_connections_from_array, color_selections,
                              select_clusters, show_cluster)
# from .interface.wordom import read_pathway
from .internal.matrix import matrix_from_interactions, matrix_from_pandas_dataframe
from .internal.procedure import draw_ciacg

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
        "-acg", nargs=1, metavar="ACGfile", help="ACG file to read (.npy)")
    parser.add_argument(
        "-rmp", nargs=1, metavar="RMPfile", help="ResidueMap file to read (.rmp)")
    arguments = parser.parse_args(argv[1:])

    # Finish pymol launch
    pymol.finish_launching(['pymol'])

    # Set variables here
    pdb = arguments.pdb[0]
    cutoffs = [float(c) for c in arguments.c]
    acg = arguments.acg[0]
    rmp = arguments.rmp[0]

    cigraph = numpy.load(acg)

    with open(rmp, 'rb') as infile:
        residuemap = pickle.load(infile)

    # Draw the loaded ciACG
    levels = draw_ciacg(cigraph, residuemap, pdb, cutoffs)


if __name__ == '__main__':
    main()
