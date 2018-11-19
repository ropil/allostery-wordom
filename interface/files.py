from pickle import dump, HIGHEST_PROTOCOL
'''
 File IO
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


def dump_pyobject(data, filename, suffix="frm"):
    """Pickles a python object, in a default manner

    :param data: data to dump
    :param filename: if not None, filename to dump into
    :param suffix: Check if suffix present in filename, otherwise add it
    """
    if filename is not None:
        outfilename = filename
        # Add proper file ending if not present
        if filename.split('.')[-1] != suffix:
            outfilename += ".{}".format(suffix)
        with open(outfilename, 'wb') as output:
            dump(data, output, HIGHEST_PROTOCOL)
