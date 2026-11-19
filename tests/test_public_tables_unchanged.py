"""
Guard the numeric values of every public residue table.

``data/public_tables_1.3.7.json`` was written by :func:`write_snapshot` from
the 1.3.7 release (master at a7abe3b), before the table parsers were
consolidated. Each table is compared through its JSON serialization, which
keeps key order, distinguishes ``4`` from ``4.0`` and round-trips floats
exactly, so any change to a value, its type or the residue order fails.
"""

import importlib
import json
from pathlib import Path

import numpy as np
import pytest

SNAPSHOT_PATH = Path(__file__).parent / "data" / "public_tables_1.3.7.json"

TABLE_MODULES = [
    "amino_acid_properties",
    "blosum",
    "pmbec",
    "residue_contact_energies",
]


def collect_public_tables():
    """
    Gather every public residue table from the table modules.

    Returns
    -------
    dict of str to object
        JSON-serializable tables keyed like ``"blosum.blosum62_dict"``.
        Dictionaries keyed by one-letter residue codes are kept as they are;
        NumPy arrays become their dtype name and nested lists of floats.
    """
    tables = {}
    for module_name in TABLE_MODULES:
        module = importlib.import_module("pepdata." + module_name)
        for name, value in sorted(vars(module).items()):
            if name.startswith("_"):
                continue
            key = "%s.%s" % (module_name, name)
            if isinstance(value, np.ndarray):
                tables[key] = {
                    "dtype": str(value.dtype),
                    "values": [[float(x) for x in row] for row in value],
                }
            elif isinstance(value, dict) and value and all(
                    isinstance(k, str) and len(k) == 1 for k in value):
                tables[key] = value
    return tables


def write_snapshot(path):
    """
    Write the tables of the installed pepdata to a JSON snapshot.

    Parameters
    ----------
    path : str or path-like
        Destination of the snapshot.
    """
    with open(path, "w") as f:
        json.dump(collect_public_tables(), f, indent=0)
        f.write("\n")


SNAPSHOT = json.loads(SNAPSHOT_PATH.read_text())


def test_snapshot_covers_the_expected_tables():
    assert len(SNAPSHOT) == 32
    modules = {key.split(".")[0] for key in SNAPSHOT}
    assert modules == set(TABLE_MODULES)


@pytest.mark.parametrize("key", sorted(SNAPSHOT))
def test_public_table_is_unchanged(key):
    current = collect_public_tables()
    assert key in current, "Public table %s was removed" % key
    assert json.dumps(current[key]) == json.dumps(SNAPSHOT[key])
