from numpy import zeros

'''
 <Decription here>
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


def matrix_from_interactions(interactions, mapping, default=0.0):
    """Generate numpy matrices from interaction dictionary

    :param interactions: dictionary read from WORDOM avgpsn
    :param mapping: residuemap, preserves residue names
    :param default: default interaction
    :return: tuple of interaction strength and frequency numpy array
    """
    # Expecting symmetric dictionary
    size = len(mapping)
    strength = zeros((size, size)) + default
    frequency = zeros((size, size)) + default
    # Populera strength och frequency via mappings från interactions
    for resa, inter in interactions.items():
        for resb, (strng, freq) in inter.items():
            a = mapping[resa] - 1
            b = mapping[resb] - 1
            strength[a][b] = strength[b][a] = strng
            frequency[a][b] = frequency[b][a] = freq
    return strength, frequency


def matrix_from_pandas_dataframe(pddframe):
    """matrix_from_pandas_dataframe returns the .values member
    
    :param pddframe: pandas dataframe
    """
    return pddframe.values
