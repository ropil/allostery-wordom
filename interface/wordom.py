import re
from collections import OrderedDict
from ..internal.map import Map
'''
 WORDOM file parsing interface
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


def read_avg_strength(infile):
    """Read residue interaction strength from PSN avg files

    :param infile: File handle pointing to WORDOM avgpsn output file
    """
    m_start = re.compile("^\*\*\* Averaged Interaction Strength \*\*\*")
    m_end = re.compile("^===")
    m_entry = re.compile("^\s*.:.\d+\s+.:.\d+\s+\d+\.\d+\s+\d+\.\d+\s*$")
    interactions = {}
    reading = False
    for line in infile:
        if reading:
            # Stop reading if end of interaction strength section
            if m_end.search(line):
                break
            else:
                if m_entry.search(line):
                    [a, b, strength, freq] = line.split()
                    if not a in interactions:
                        interactions[a] = {}
                    if not b in interactions:
                        interactions[b] = {}
                    # Assign symmetrically
                    interactions[a][b] = interactions[b][a] = (float(strength), float(freq))
        # Start reading when header found
        elif m_start.search(line):
            reading = True
    return interactions


def read_avg_clusters(infile):
    """ Read clusters from PSN avg files

    :param infile: File handle pointing to WORDOM avgpsn output file
    """
    m_start = re.compile("^\*\*\* Stable Cluster Compositions \*\*\*")
    m_imin = re.compile("^Imin:")
    m_freq = re.compile("^Freq:")
    m_end = re.compile("^===")
    m_entry = re.compile("^C\s*\d+:")
    reading = False
    clusters = {}
    current_imin = None
    current_freq = None
    for line in infile:
        if reading:
            if m_end.search(line):
                return clusters
            else:
                if m_entry.search(line):
                    entries = ":".join(line.split(':')[1:])
                    clusters[current_imin][current_freq].append(entries.split())
                elif m_imin.search(line):
                    current_imin = float(line.split()[1])
                    clusters[current_imin] = {}
                elif m_freq.search(line):
                    current_freq = float(line.split()[1])
                    clusters[current_imin][current_freq] = []
        elif m_start.search(line):
            reading = True
    return clusters


def read_avg_residuemap(infile):
    """ Read sequence definition from PSN avg file, returning sequence Map

    :param infile: File handle pointing to WORDOM avgpsn output file
    :return: Returns an internal.map.Map object mapping the .pdb
             residues to WORDOM id's from "Seq" section of the avgpsn-file
    """
    m_start = re.compile("^\*\*\* Seq \*\*\*")
    m_end = re.compile("^============")
    m_entry = re.compile("^\s*\d+\s+.:\d+\s+.:.\d+\s*$")
    residuemap = OrderedDict()
    reading = False
    for line in infile:
        if reading:
            # Stop reading if end of interaction strength section
            if m_end.search(line):
                break
            else:
                if m_entry.search(line):
                    [num, resname, reslongname] = line.split()
                    residuemap[resname] = int(num) - 1
        # Start reading when header found
        elif m_start.search(line):
            reading = True
    return residuemap


def get_chain_offsets(chainlist, chainlength, chainpadding):
    """Generates chain offsets

    :param chainlist: Sorted list of chain letters
    :param chainlength: chain length dictionary
    :param chainpadding: chain padding dictionary
    :return: dictionary of chain offsets
    """
    chainoffset = {}
    for a in range(len(chainlist)):
        if a > 0:
            chainoffset[chainlist[a]] = sum([0] + [chainlength[chainlist[b]] +
                                                   chainpadding[chainlist[b]]
                                                   for b in range(a)])
        else:
            chainoffset[chainlist[a]] = 0
    return chainoffset


def wordom_to_map(residuemap, chainlength=None):
    """Generates a WORDOM sequence to contact map sequence


    :param residuemap: OrderedDict, mapping entries of type "A:21" to int
    :param chainlength: dictionary with chain letters as keys and int lengths
                        as items
    :return: tuple of; Map for WORDOM sequence to contact map numbering, total
             contact map size, chainpadding dict, chainoffset dict, chainlength
             dict.
    """
    chainstart = {}
    # Get chain sizes, maximum and minimum PDB residue number
    if chainlength is None:
        chainlength = {}
        for entry in residuemap:
            [chain, num] = entry.split(':')
            # Initialize empty entries
            if not chain in chainlength:
                chainlength[chain] = int(num)
                chainstart[chain] = int(num)
            # get maximum and minimum
            chainlength[chain] = max([int(num), chainlength[chain]])
            chainstart[chain] = min([int(num), chainlength[chain]])
    # Calculate chain padding for negative PDB indices
    chainpadding = {}
    for chain in chainstart:
        if chainstart[chain] < 1:
            chainpadding[chain] = 1 - chainstart[chain]
        else:
            chainpadding[chain] = 0
    # Get an alphabetically sorted list of chains
    chainlist = list(chainlength.keys())
    chainlist.sort()
    # Calculate the chain offsets, padding for negative starting indices
    chainoffset = get_chain_offsets(chainlist, chainlength, chainpadding)
    # chainoffset = {}
    # for a in range(len(chainlist)):
    #     if a > 0:
    #         chainoffset[chainlist[a]] = sum([0] + [chainlength[chainlist[b]] +
    #                                                chainpadding[chainlist[b]]
    #                                                for b in range(a)])
    #     else:
    #         chainoffset[chainlist[a]] = 0
    # Generate the residuemap to map mapping
    cmapseq = []
    for entry in residuemap:
        [chain, num] = entry.split(':')
        # Consider the chain padding and offset; numbering 1 in PDB as 0
        cmapseq.append(chainpadding[chain] + chainoffset[chain] + int(num) - 1)
    # Return the map, padding and contact map size
    cmapsize = sum([chainpadding[i] for i in chainpadding]) + sum(
        [chainlength[j] for j in chainlength])
    return Map(cmapseq), cmapsize, chainpadding, chainoffset, chainlength

