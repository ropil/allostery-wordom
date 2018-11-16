#!/usr/bin/env python3
# import pymol
import numpy
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
            # Expect residues on format C:AX, with C & A char and X int
            residues.append("(chain {} and resi {})".format(
                resi.split(':')[0], resi.split(':')[1][1:]))
        residues = " or ".join(residues)
        clusternames.append("c{:>02d}".format(cnum))
        selection = "({}) and name CA".format(residues)
        cmd.select(clusternames[-1], selection)
    return clusternames


def bond_connections(clusters, interactions):
    # Quickfix - find min and max
    # this filtering should be done before supplying to this function
    # That is, the "interactions" data should be already formatted
    # as a float in a way that makes it look nice with bond 
    # stick_radius
    minimum = None
    maximum = None
    for cluster in clusters:
        for i in range(len(cluster)):
            for j in range(i, len(cluster)):
                resa = cluster[i]
                resb = cluster[j]
                if resa in interactions:
                    if resb in interactions[resa]:
                        strength = interactions[resa][resb][0]
                        if minimum is None:
                            minimum = strength
                        else:
                            minimum = strength if minimum > strength else minimum
                        if maximum is None:
                            maximum = strength
                        else: 
                            maximum = strength if maximum < strength else maximum
    for cluster in clusters:
        for i in range(len(cluster)):
            for j in range(i, len(cluster)):
                resa = cluster[i]
                resb = cluster[j]
                # Only draw bonds if interaction strength at all present
                if resa in interactions:
                    if resb in interactions[resa]:
                        a = "chain {} and resi {} and name CA".format(
                            resa.split(':')[0], resa.split(':')[1][1:])
                        b = "chain {} and resi {} and name CA".format(
                            resb.split(':')[0], resb.split(':')[1][1:])
                        cmd.bond(a, b)
                        strength = 1.0 if maximum == minimum else 0.1 + (0.9 * ((interactions[resa][resb][0] - minimum) / (maximum - minimum)))
                        cmd.set_bond("stick_radius", strength, a, b)


def bond_connections_from_array(interactiongraph, residuemap, cutoff=0.0):
    # Make strength be proportional to girth of bonds
    graph = numpy.sqrt(numpy.absolute(interactiongraph))
    # Simple min-max scaling parameters
    minimum = numpy.min(graph)
    maximum = numpy.max(graph)
    residues = list(residuemap.keys())
    shown = []
    for i in range(graph.shape[0]):
        for j in range(i, graph.shape[1]):
            # Only draw bonds if interaction strength over cutoff threshold
            if graph[i][j] > cutoff:
                resa = residues[i]
                resb = residues[j]
                shown.append(resa)
                shown.append(resb)
                a = "chain {} and resi {} and name CA".format(
                    resa.split(':')[0], resa.split(':')[1][1:])
                b = "chain {} and resi {} and name CA".format(
                    resb.split(':')[0], resb.split(':')[1][1:])
                cmd.bond(a, b)
                strength = 1.0 if maximum == minimum else 0.1 + (0.9 * ((graph[i][j] - minimum) / (maximum - minimum)))
                cmd.set_bond("stick_radius", strength, a, b)
    return shown


def bond_colors_from_array(colorarray, residuemap, cutoff=0.0, colorprefix="path_"):
    residues = list(residuemap.keys())
    colored = []
    colors = {}
    colorindex = 0
    # Expect rgb channels over first dimension
    for i in range(colorarray.shape[1]):
        for j in range(i, colorarray.shape[2]):
            # Only draw bonds if interaction strength over cutoff threshold
            resa = residues[i]
            resb = residues[j]
            colored.append((resa, resb))
            a = "chain {} and resi {} and name CA".format(
                resa.split(':')[0], resa.split(':')[1][1:])
            b = "chain {} and resi {} and name CA".format(
                resb.split(':')[0], resb.split(':')[1][1:])
            #cmd.bond(a, b)
            color = tuple(colorarray[0:3,i,j])
            if color not in colors:
                colorindex += 1
                colors[color] = colorindex
                cmd.set_color("{}{}".format(colorprefix, colorindex), color)
            else:
                colorindex = colors[color]
            colorname = "{}{}".format(colorprefix, colorindex)
            print("Applying color {}, named as {}, to residues {}".format(color, colorname, colored[-1]))
            cmd.set_bond("stick_color", colorname, a, b)
    return colored, colors


def show_cluster(clusters):
    for cluster in clusters:
        for resi in cluster:
            sele = "chain {} and resi {} and name CA".format(
                resi.split(':')[0], resi.split(':')[1][1:])
            cmd.show("sticks", sele)
            cmd.show("spheres", sele)

