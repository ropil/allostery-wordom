#!/usr/bin/env python3
if __name__ == "__main__" and __package__ is None:
    __package__ = "allostery-wordom"

import pymol
from pymol import cmd

from .interface.pymol import (bond_connections_from_array, color_selections,
                              select_clusters, show_cluster)
from .interface.wordom import read_avg_strength, read_avg_residuemap, read_correlations
from .internal.map import Map
from .internal.matrix import matrix_from_interactions, matrix_from_pandas_dataframe

import numpy
import matplotlib.pyplot as plt
from pandas import DataFrame
'''
 Display the icACG in an interactive PyMOL session
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
        "allostery_icacg_pymol  Copyright (C) 2018  Robert Pilstål;",
        "This program comes with ABSOLUTELY NO WARRANTY.",
        "This is free software, and you are welcome to redistribute it",
        "under certain conditions; See the supplied Apache License,",
        "Version 2.0 for the specific language governing permissions",
        "and limitations under the License.",
        "    http://www.apache.org/licenses/LICENSE-2.0"
    ])


# Main; for callable scripts
def main():
    from argparse import ArgumentParser
    from sys import argv
    parser = ArgumentParser(
        description= "Open and display a WORDOM PSN and cross correlation " +
                     "analysis as an interaction-correlation Allosteric " +
                     "Communication Graph (icACG) in PyMOL."
    )
    parser.add_argument(
        "-i",
        nargs=1,
        default=[None],
        metavar="float",
        help="Imin to use (must be present in AVGfile)" +
        ", default=Use lowest found")
    parser.add_argument(
        "-f",
        nargs=1,
        default=[None],
        metavar="float",
        help="Freq to use (must be present for selected Imin)" +
        ", default=Use lowest found")
    parser.add_argument(
        "-c",
        nargs=1,
        default=[0.0],
        metavar="float",
        help="Cutoff to use for showing connections in graph" +
        ", default=0.0")
    parser.add_argument(
        "-avg", nargs=1, metavar="AVGfile", help="Wordom PSN avg file")
    parser.add_argument(
        "-cor",
        nargs=1,
        metavar="CORRfile",
        help="Wordom cross-correlation file")
    parser.add_argument(
        "-pdb", nargs=1, metavar="PDBfile", help="PDB file to draw")
    parser.add_argument("-plot", action="store_true", default=False, help="Plot icACG value distribution")
    arguments = parser.parse_args(argv[1:])

    # Finish pymol launch
    pymol.finish_launching(['pymol'])

    # Set variables here
    pdb = arguments.pdb[0]
    avg = arguments.avg[0]
    cor = arguments.cor[0]
    cutoff = float(arguments.c[0])
    imin = arguments.i[0]
    freq = arguments.f[0]
    icplot = arguments.plot

    interactions = {}
    with open(avg, 'r') as infile:
        residuemap = read_avg_residuemap(infile)
        mapping = Map([int(i.split(':')[-1][1:]) for i in residuemap.keys()], inverse_sequence=list(range(len(residuemap))))
        infile.seek(0)
        interactions = read_avg_strength(infile)

    (strength, frequency) = matrix_from_interactions(interactions, residuemap)

    with open(cor, 'r') as infile:
        correlation_table = read_correlations(infile)

    correlation = matrix_from_pandas_dataframe(correlation_table)

    icgraph = strength * correlation

    if icplot:
         plt.figure()
         df = DataFrame({'a': icgraph.reshape(icgraph.shape[0] * icgraph.shape[1])}, columns=['a'])
         df.plot.hist(stacked=True)
         plt.show()

    cmd.load(pdb)
    cmd.hide("everything")
    cmd.show("ribbon")

    # Create bindings and selections, and color them
    residues = bond_connections_from_array(icgraph, residuemap, cutoff=cutoff)
    selections = select_clusters([residues])
    colors = color_selections(selections)

    # Show clusters
    show_cluster([residues])


if __name__ == '__main__':
    main()
