#!/usr/bin/env python3
import pymol
from colorsys import hsv_to_rgb
from pymol import cmd
'''
 PyMOL interface, accessing API but not redistributing PyMOL source
 Allostery-WORDOM PyMOL interface Copyright (C) 2015-2018  Robert Pilstål
 Open-Source PyMOL is Copyright (C) Schrodinger, LLC.

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


# Library functions


# Colors a range of selections (need to be created first)
def color_selections(selections):
    cl = len(selections)
    hsv_range = [(float(i + 1) * 1.0 / float(cl), 0.5, 0.5)
                 for i in range(cl)]
    # *x means the same as in
    # rgb_range = map(lambda x: hsv_to_rgb(x[0], x[1], x[2]), hsv_range)
    rgb_range = map(lambda x: hsv_to_rgb(*x), hsv_range)
    colornames = []
    for [selection, color] in zip(selections, rgb_range):
        cmd.set_color(selection, color)
        colornames.append(selection)
        cmd.color(colornames[-1], selection)
    return colornames
cmd.extend("color_selections", color_selections)


def select_clusters(clusters):
    clusternames = []
    for [cluster, cnum] in zip(clusters, range(len(clusters))):
        residues = []
        for resi in cluster:
            residues.append("(chain {} and resi {})".format(
                resi.split(':')[0], resi.split(':')[1]))
        residues = " or ".join(residues)
        clusternames.append("c{:>02d}".format(cnum))
        cmd.select(clusternames[-1], "({}) and name CA".format(residues))
    return clusternames


def bond_connections(clusters, interactions):
    for cluster in clusters:
        for i in range(len(cluster)):
            for j in range(i, len(cluster)):
                resa = cluster[i]
                resb = cluster[j]
                # Only draw bonds if interaction strength at all present
                if resa in interactions:
                    if resb in interactions[resa]:
                        a = "chain {} and resi {} and name CA".format(
                            resa.split(':')[0], resa.split(':')[1])
                        b = "chain {} and resi {} and name CA".format(
                            resb.split(':')[0], resb.split(':')[1])
                        cmd.bond(a, b)
                        cmd.set_bond("line_width",
                                     1.0 + interactions[resa][resb],
                                     a, b)


def show_cluster(clusters):
    for cluster in clusters:
        for resi in cluster:
            sele = "chain {} and resi {} and name CA".format(
                resi.split(':')[0], resi.split(':')[1])
            cmd.show("lines", sele)
            cmd.show("spheres", sele)

