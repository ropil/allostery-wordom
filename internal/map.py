from copy import copy
from collections import OrderedDict
import sqlite3
'''
 OrderedDict with extra features for protein seq-seq mappings
 Copyright (C) 2016-2018  Robert Pilstål

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


class Map(OrderedDict):
    def __init__(self,
                 sequence,
                 *args,
                 len=None,
                 inverse=None,
                 inverse_sequence=None,
                 inverse_len=None,
                 **kwargs):
        """Maps a sequence of integers to another sequence of integers, based
        on the OrderedDict class. Gaps are coded as mappings to None

        :param sequence: an iterable sequence of integers, 0 will map t
                         sequence[0], 1 will map to sequence[1] and so on. Gaps
                         encoded as None.
        :param args:
        :param inverse: Specify the inverse of the map, otherwise make it from
                        the provided sequence.
        :param kwargs:
        """
        OrderedDict.__init__(self, *args, **kwargs)
        for t in enumerate(sequence):
            self[t[0]] = t[1]
        # Initiate the inverse
        if inverse is None:
            if inverse_sequence is None:
                # We have gotten no idea, so we guess
                self.inverse = Map(self.guess_inverse(sequence), inverse=self)
            else:
                # We have gotten the inverse sequence,
                # so initiate a map with it
                self.inverse = Map(inverse_sequence, inverse=self)
        else:
            # We have gotten an map, so we just apply it
            self.inverse = inverse

    def __conform__(self, protocol):
        if protocol is sqlite3.PrepareProtocol:
            return "\n".join([str(self), str(self.inverse)])

    def __str__(self):
        return self.draw_sequential_alignment()

    def push(self, a, b):
        # Map maps a to b, get a copy of b with a substituted where it maps
        # Make a new copy of b
        c = [copy(b[i]) for i in range(len(b))]
        for t in enumerate(a):
            # If index maps to
            if t[0] in self:
                c[t[0]] = copy(a[self[t[0]]])
        return c

    def fetch(self, a, b):
        # Map maps a to b, get a copy of a with b substituted where it maps
        # Make a new copy of b
        c = []
        for t in enumerate(a):
            # If index maps to
            if t[0] in self:
                c.append(copy(b[self[t[0]]]))
            else:
                c.append(copy(t[1]))
        return c

    def gapped(self, a, size, gap='-'):
        # Create a gapped sequence
        b = [copy(gap) for i in range(size)]
        return self.push(a, b)

    def gapslice(self, a, gap='-'):
        # Get slice size
        imax = max(self.items())
        imin = min(self.items())
        # Create a gapped sequence slice
        c = [copy(gap) for i in range(imax - imin)]
        # Insert a on mapped positions in c
        for t in enumerate(a):
            if t[0] in self:
                c[self[t[0]] - imin] = copy(t[1])
        return c

    def guess_inverse(self, sequence):
        # Guess the inverse from sequence
        # Create template of same length as sequence
        reordered = [None] * max(sequence)
        for pos in enumerate(sequence):
            if pos[1] is not None:
                reordered[pos[1]] = pos[0]
        return reordered

    def domain(self, a):
        c = []
        # Return an ordered copy of the mapped elements of a
        for i in self:
            c.append(copy(a[self[i]]))
        return c

    def draw_sequential_alignment(self):
        # This is assigned backwards and will only work if alignments are
        # both sequential.
        pos = 0
        seq = ""
        for i in self.inverse:
            current = self.inverse[i]
            # if a residue aligns
            if current is not None:
                # if self.inverse[i] > pos:
                #     # Insert any positionals that are inserts pertaining to
                #     # other
                #     seq += "X" * (i - pos)
                # else:
                #     seq += "X"
                seq += "X" * (current - pos + 1)
                pos = current + 1
            else:
                seq += "-"
        # Append any trailing parts of this sequence
        if pos < len(self) - 1:
            seq += "X" * (len(self) - pos)
        return seq

    def picture(self, a):
        # Return the domain picture of a
        """Get the picture of domain a, sorted in the mapped space

        :param a: sequence to act on
        :return: domain of a, sorted by the mapped order
        """
        return [copy(a[self.inverse[i]]) for i in self.inverse]
