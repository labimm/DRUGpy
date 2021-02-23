import hashlib
import random
import shlex
import subprocess
import tempfile
from contextlib import contextmanager
import os

import numpy as np
import pandas as pd
import scipy.spatial
import seaborn as sb
from matplotlib import pyplot as plt
from pymol import cmd as pm
from pymol import stored
from pymol.plugins import pref_get


PLUGIN_DATA_DIR = os.path.expanduser("~/.pymol/labimm")
os.makedirs(PLUGIN_DATA_DIR, exist_ok=True)

#
# Generic
#


def pairwise(iter1, iter2):
    """Pairwise element iteration."""
    for item1 in iter1:
        for item2 in iter2:
            yield (item1, item2)


#
# PyMOL specific
#

@contextmanager
def disable_feedback(what, level):
    """Disable feedback."""
    pm.feedback("disable", what, level)
    yield
    pm.feedback("enable", what, level)


@contextmanager
def settings(**kwargs):
    """Set options temporarily."""
    orig = {}
    for key, value in kwargs.items():
        orig[key] = pm.get(key)
        pm.set(key, value)
    yield
    for key, value in orig.items():
        pm.set(key, value)


#
# Atomic helpers
#


def contact_matrix(a, b, radius=None):
    """
    Compute the contact matrix between A and B.
    When radius is None the contact is defined by:
        Distance(a, b) < VDW(a) + VDW(b)
    When radius is a number the contact is defined by:
        Distance(a, b) < radius
    """
    d = scipy.spatial.distance_matrix(a[["x", "y", "z"]], b[["x", "y", "z"]])
    if radius is None:
        return ((d - b["vdw"].values).T - a["vdw"].values).T <= 0
    else:
        return (d - radius) <= 0


def get_atoms(sel, attrs, state=1):
    """Get the atoms and attributes of a selection."""
    coords = None
    if "coords" in attrs:
        coords = pm.get_coords(sel, state)
        attrs.remove("coords")
    atoms = pd.DataFrame(coords, columns=["x", "y", "z"])

    if attrs:
        fields_str = ", ".join(attrs)
        stored.atoms = []
        pm.iterate_state(state, sel, f"stored.atoms.append(({fields_str}))")
        atoms = pd.concat([atoms, pd.DataFrame(stored.atoms, columns=attrs)], axis=1)
        del stored.atoms
    return atoms


def count_molecules(selection="all"):
    """
    By Thomas Holder.
    """
    tmpsele = pm.get_unused_name("_tmp")
    count = 0
    if pm.select(tmpsele, selection):
        count += 1
        while pm.select(tmpsele, f"{tmpsele} &! bm. first {tmpsele}"):
            count += 1
    pm.delete(tmpsele)
    return count


#
# Etcetera
#


def run(command, **kwargs):
    if isinstance(command, str):
        command = shlex.split(command)
    ret = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kwargs)
    output = ret.stdout.decode()
    success = ret.returncode == 0
    return output, success


def rscript(script):
    tmp_fd, tmp_path = tempfile.mkstemp(".R")
    with open(tmp_fd, "w") as tmp_file:
        tmp_file.write(script)
    return run([pref_get("LABIMM_RSCRIPT"), tmp_path])


#
# Commands
#


def fo_(sel1, sel2, radius=2, state1=1, state2=1):
    atoms1 = get_atoms(sel1, ["coords", "elem", "vdw"], state1)
    atoms2 = get_atoms(sel2, ["coords", "elem", "vdw"], state2)

    a = atoms1[atoms1["elem"] != "H"]
    b = atoms2[atoms2["elem"] != "H"]
    contacts = np.any(contact_matrix(a, b, radius), axis=1)
    num_contacts = np.sum(contacts)
    total_atoms = len(a)
    fo = num_contacts / total_atoms
    print(f"  Fractional Overlap = {fo}")
    return fo


@pm.extend
def fo(sel1, sel2, radius=2, state1=1, state2=1):
    """
    Compute the fractional overlap of sel1 respective to sel2.
        FO = Nc/Nt

    Nc is the number of atoms of sel1 in contact with sel2. Nt is the number of atoms
    of sel1.

    Hydrogen atoms are ignored.

    If the contact radius is 0 then the VdW radii will be used.

    The states are for select a single state from a multi-state objects.

    OPTIONS:
        sel1    ligand object.
        sel2    hotspot object.
        radius  the radius so two atoms are in contact (default: 2).
        state1  state of sel1.
        state2  state of sel2.

    EXAMPLES:
        fo ref_lig, ftmap1234.D.003
        fo ref_lig, ftmap1234.CS.000_016
    """
    print(sel1 + " / " + sel2)
    fo_(sel1, sel2, radius, state1, state2)
    print(sel2 + " / " + sel1)
    fo_(sel2, sel1, radius, state2, state1)
