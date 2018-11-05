#!/usr/bin/env python3
import pymol
from pymol import cmd

from .interface.pymol import (bond_connections, color_selections,
                              select_clusters, show_cluster)
from .interface.wordom import read_avg_clusters, read_avg_strength


'''
 Display the avgpsn analysis in an interactive PyMOL session
 Copyright (C) 2015-2018  Robert Pilstål

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
        "allostery_psn_pymol  Copyright (C) 2015-2018  Robert Pilstål;",
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
        description="Open and display a WORDOM PSN cluster analysis in PyMOL.")
    parser.add_argument(
        "-c", nargs=1, default=[None], metavar="float",
        help="Cutoff to use (must be present in AVGfile)" +
        ", default=Use lowest found")
    parser.add_argument(
        "-show", nargs=1, default=[None], metavar="int[,int[...]]",
        help="Show specified clusters")
    parser.add_argument(
        "avg", nargs=1, metavar="AVGfile",
        help="Wordom PSN avg file")
    parser.add_argument(
        "pdb", nargs=1, metavar="PDBfile", help="PDB file to draw")
    arguments = parser.parse_args(argv[1:])

    # Finish pymol launch
    pymol.finish_launching()

    # Set variables here
    pdb = arguments.pdb[0]
    avg = arguments.avg[0]
    cut = arguments.c[0]
    shw = arguments.show[0]

    interactions = {}
    clusters = {}
    with open(avg, 'r') as infile:
        interactions = read_avg_strength(infile)
        infile.seek(0)
        clusters = read_avg_clusters(infile)

    # Select the cutoff
    if cut is not None:
        # If provided
        cut = float(cut)
    else:
        # Default
        cutoffs = list(clusters.keys())
        cutoffs.sort()
        cut = cutoffs[0]

    # Select clusters
    clusters = clusters[cut]

    cmd.load(pdb)
    cmd.hide("everything")
    cmd.show("ribbon")

    # Create bindings and selections, and color them
    bond_connections(clusters, interactions)
    selections = select_clusters(clusters)
    colors = color_selections(selections)

    # Show clusters
    if shw is None:
        show_cluster(clusters)
    else:
        shw = [int(c) for c in shw.split(',')]
        for c in shw:
            show_cluster(clusters[c - 1])


if __name__ == '__main__':
    main()
