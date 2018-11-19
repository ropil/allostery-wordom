from numpy import add, around, divide, multiply, ones, subtract, zeros
from pandas import DataFrame

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


def align_dataframes(framea, frameb, fill_value = 0.0):
    """Use pandas DataFrame structure to align two-dimensional data

    :param framea: First pandas dataframe to align
    :param frameb: Other pandas dataframe to align
    :param fill_value: default fill value (0.0 float)
    return: tuple of aligned frames
    """
    zeroframe = frameb.copy()
    zeroframe[:] = fill_value
    aligneda = framea.add(zeroframe, fill_value = fill_value)
    zeroframe = framea.copy()
    zeroframe[:] = fill_value
    alignedb = frameb.add(zeroframe, fill_value = fill_value)
    return aligneda, alignedb


def dataframe_from_dictionary(data, indexmap = None, fillna = None):
    """Create a pandas dataframe from dictionary of dicts

    :param data: dictionary of dicts
    :param indexmap: a dictionary for mapping indices to new names,
                     if None, no mapping will be done (default)
    :param fillna: value to fill missing values with. If None, no
                   filling will be done (default)
    :return: pandas dataframe
    """
    # Convert dictionary to DataFrame
    df = DataFrame.from_dict(data,orient='index')
    if indexmap is not None:
        df = df.rename(index = indexmap, columns = indexmap)
    if fillna is not None:
        df = df.fillna(value = fillna)
    return df


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


def matrix_to_colorarray(matrix, hue_from = (1.0, 1.0, 1.0), hue_to = (1.0, 0.0, 0.0), channel_max = 255.0):
    """Converts a matrix into a color channel array

    :param matrix: matrix to convert
    :param hue_from: tuple of hue-values, starting color for null
                     strength. Default three with first channel
                    1.0 and the rest 0.
    :param hue_to: tuple with ending hue-values, for max strength.
                   Coloring will be interpolated between hue_from
                   and hue_to.
    :param channel_max: Maximum level of channel, default 255.0
    :return: array with first dimension being number of channels,
             providing one matrix per channel over the other
             dimensions
    """
    rgb_matrix = zeros((len(hue_to), matrix.shape[0], matrix.shape[1]))

    # Find max, min and span
    element_min = matrix.min()
    element_max = matrix.max()
    element_range = element_max - element_min

    print("Minimum = {}, Maximum = {}, span = {}".format(element_min, element_max, element_range))

    # Normalize matrix between 0.0 and 1.0
    normed_matrix = divide(subtract(matrix, element_min), element_range)

    # Calculate rgb channels using specified "hue", rounding off
    for i in range(len(hue_to)):
        rgb_matrix[i,:,:] = around(multiply(add(multiply(normed_matrix, hue_to[i] - hue_from[i]), multiply(ones((matrix.shape[0], matrix.shape[1])), hue_from[i])), channel_max))
        #rgb_matrix[i,:,:] = around(multiply(normed_matrix, hue[i] * channel_max))

    return rgb_matrix
