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

from functools import wraps

def _prepare_memoization_key(args, kwargs):
    """
    Make a tuple of arguments which can be used as a key
    for a memoized function's lookup_table. If some object can't be hashed
    then used its __repr__ instead.
    """
    key_list = []
    for arg in args:
        try:
            hash(arg)
            key_list.append(arg)
        except:
            key_list.append(repr(arg))
    for (k, v) in kwargs.items():
        try:
            hash(k)
            hash(v)
            key_list.append((k, v))
        except:
            key_list.append((repr(k), repr(v)))
    return tuple(key_list)

def memoize(fn):
    lookup_table = {}

    @wraps(fn)
    def wrapped_fn(*args, **kwargs):
        key = _prepare_memoization_key(args, kwargs)
        if key not in lookup_table:
            lookup_table[key] = fn(*args, **kwargs)
        return lookup_table[key]

    return wrapped_fn
