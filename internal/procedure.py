from pymol import cmd
from ..interface.pymol import bond_connections_from_array, select_clusters, color_selections, show_cluster

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
