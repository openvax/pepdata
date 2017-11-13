# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from __future__ import print_function, division, absolute_import

class AminoAcid(object):
    def __init__(
            self, full_name, short_name, letter, contains=None):
        self.letter = letter
        self.full_name = full_name
        self.short_name = short_name
        if not contains:
            contains = [letter]
        self.contains = contains

    def __str__(self):
        return (
            ("AminoAcid(full_name='%s', short_name='%s', letter='%s', "
             "contains=%s)") % (
            self.letter, self.full_name, self.short_name, self.contains))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return other.__class__ is AminoAcid and self.letter == other.letter
