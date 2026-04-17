"""
Create and Analyze System Stacking Fingerprints (SSFs)

This module contains the functions to create and analyze System
Stacking Fingerprints. It allows the user to generate pairwise
distance data across a trajectory and compare these matrices across
different trajectories.
"""

import mdtraj as md
import numpy as np
from numpy import typing
from .residue_movement import calc_center_3pts
from .vector import *
from .file_manipulation import SmartIndexingAction
from .visualization import NoResidues, create_axis_labels, display_arrays_as_video
import sys
import concurrent.futures
import functools
import math
from numba import njit, prange
import multiprocessing as mp
from collections.abc import Mapping

_NUCLEOTIDE_NAMES = {"A", "A5", "A3", "G", "G5", "G3", "C", "C5", "C3",
                     "T", "T5", "T3", "U", "U5", "U3", "INO"}
_PYRIMIDINE_NAMES = {"C", "C5", "C3","T" "T5", "T3", "U", "U5", "U3"}
_PURINE_INO_NAMES = {"A", "A5", "A3", "G", "G5", "G3","INO"}
_AROMATIC_AA_NAMES = {"PHE", "TYR", "TRP", "HIS", "HID", "HIE", "HIP"}
_ARG_AA_NAMES = {"ARG"}
_PHE_NAMES = {"PHE"}
_TYR_NAMES = {"TYR"}
_TRP_NAMES = {"TRP"}
_HIS_NAMES = {"HIS","HIE","HID","HIP"}
_ARG_NAMES = {"ARG"}

_PYRIMIDINE_RING_ATOMS = {"N1","C2","N3","C4","C5","C6"}
_PURINE_INO_RING_ATOMS = {"N1","C2","N3","C4","C5","C6","N7","C8","N9"}
_PHE_RING_ATOMS = {"CG","CD1","CE1","CZ","CE2","CD2"}
_TYR_RING_ATOMS = {"CG","CD1","CE1","CZ","CE2","CD2"}
_TRP_RING_ATOMS = {"CG","CD1","NE1","CE2","CD2","CE3","CZ3","CH2","CZ2"}
_HIS_RING_ATOMS = {"CG","ND1","CE1","NE2","CD2"}
_ARG_ATOMS = {"NE","CZ","NH1","NH2"}


def calculate_residue_distance(trj: md.Trajectory, 
                               res1: int, 
                               res2: int, 
                               res1_atoms: tuple = ("C2", "C4", "C6"),
                               res2_atoms: tuple = ("C2", "C4", "C6"),
                               frame: int = 1) -> Vector:
    """
    Calculates the vector between two residues with x, y, z units in Angstroms.

    Calculates the distance between the center of two residues. The center is defined
    by the average x, y, z position of three passed atoms for each residue (typically
    every other carbon on the 6-carbon ring of the nucleotide base).

    Parameters
    ----------
    trj : md.Trajectory
        Single frame trajectory.
    res1 : int
        1-indexed residue number of the first residue (PDB Column 5).
    res2 : int
        1-indexed residue number of the second residue (PDB Column 5).
    res1_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 1.
    res2_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 2.
    frame : int, default=1
        1-indexed frame number of trajectory to calculate the distance.

    Returns
    -------
    distance_res12 : Vector
        Vector from the center of geometry of residue 1 to residue 2.
    
    See Also
    --------
    get_residue_distance_for_frame : Calculates pairwise distances between all residues in a given frame.

    Examples
    --------
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                              'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                              residues = {426,427}, 
    ...                              atoms = {'C2','C4','C6'})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    >>> distance_vec = st.calculate_residue_distance(
    ...     trajectory=filtered_traj, 
    ...     res1_num=426, 
    ...     res2_num=427, 
    ...     res1_atoms=("C2", "C4", "C6"), 
    ...     res2_atoms=("C2", "C4", "C6"), 
    ...     frame=1
    ... )
    >>> distance_vec.magnitude()
    7.5253396
    """
    trj = trj[frame-1]

    # Correct for mdtraj 0-indexing
    res1 = res1 - 1 
    res2 = res2 - 1

    topology = trj.topology
    res1_atom_indices = topology.select("resSeq " + str(res1))
    res2_atom_indices = topology.select("resSeq " + str(res2))
    res1_name = topology.atom(res1_atom_indices[0]).residue.name
    res2_name = topology.atom(res2_atom_indices[0]).residue.name

    if (res1_name not in _NUCLEOTIDE_NAMES) or (res2_name not in _NUCLEOTIDE_NAMES):
        return Vector(0,0,0)
    
    desired_res1_atom_indices = topology.select("(name " + res1_atoms[0] + " or name " + res1_atoms[1] + " or name " + res1_atoms[2] + ") and residue " + str(res1))
    desired_res2_atom_indices = topology.select("(name " + res2_atoms[0] + " or name " + res2_atoms[1] + " or name " + res2_atoms[2] + ") and residue " + str(res2))

    # convert nanometer units in trajectory.xyz to Angstroms
    res1_atom_xyz = trj.xyz[0, desired_res1_atom_indices, :] * 10
    res2_atom_xyz = trj.xyz[0, desired_res2_atom_indices, :] * 10
    vectorized_res1_atom_xyz = [Vector(x,y,z) for [x,y,z] in res1_atom_xyz]
    vectorized_res2_atom_xyz = [Vector(x,y,z) for [x,y,z] in res2_atom_xyz]
    res1_center_of_geometry = calc_center_3pts(*vectorized_res1_atom_xyz)
    res2_center_of_geometry = calc_center_3pts(*vectorized_res2_atom_xyz)

    distance_res12 = res2_center_of_geometry - res1_center_of_geometry
    return distance_res12

def determine_default_ring_atoms(res_name: str) -> set:
    if res_name in _PYRIMIDINE_NAMES:
        return _PYRIMIDINE_RING_ATOMS
    elif res_name in _PURINE_INO_NAMES:
        return _PURINE_INO_RING_ATOMS
    elif res_name in _PHE_NAMES:
        return _PHE_RING_ATOMS
    elif res_name in _TYR_NAMES:
        return _TYR_RING_ATOMS
    elif res_name in _TRP_NAMES:
        return _TRP_RING_ATOMS
    elif res_name in _HIS_NAMES:
        return _HIS_RING_ATOMS
    elif res_name in _ARG_NAMES:
        return _ARG_ATOMS
    else:
        return set()

# --- core helpers used by BOTH frame() and trajectory() ---
@njit(parallel=True, fastmath=True)
def _pairwise_dist_numba(centers):
    n = centers.shape[0]
    out = np.zeros((n, n), dtype=np.float32)
    for i in prange(n):
        xi0, xi1, xi2 = centers[i, 0], centers[i, 1], centers[i, 2]
        if np.isnan(xi0):
            continue
        for j in range(i + 1, n):
            xj0, xj1, xj2 = centers[j, 0], centers[j, 1], centers[j, 2]
            if np.isnan(xj0):
                continue
            dx = xj0 - xi0; dy = xj1 - xi1; dz = xj2 - xi2
            d = (dx*dx + dy*dy + dz*dz) ** 0.5
            out[i, j] = d; out[j, i] = d
    return out

def _build_ring_index_map(topology, atoms_spec=None):
    """
    atoms_spec can be:
      - None: use determine_default_ring_atoms(res.name) for each residue
      - set[str]: same atom set for all residues (legacy behavior)
      - dict[str, set[str]]: per-residue atom selections
    """
    idx_map = []

    for res in topology.residues:
        # Case 1: per-residue mapping
        if isinstance(atoms_spec, Mapping):
            names = atoms_spec.get(res.name)
            if names is None:
                names = determine_default_ring_atoms(res.name)

        # Case 2: single global atom set (legacy)
        elif isinstance(atoms_spec, set):
            names = atoms_spec or determine_default_ring_atoms(res.name)

        # Case 3: no override -> pure default behavior
        elif atoms_spec is None:
            names = determine_default_ring_atoms(res.name)

        else:
            raise TypeError("atoms_spec must be None, a set, or a mapping of resName -> set[atomNames]")

        # Build atom index array
        idx = np.array([a.index for a in res.atoms if a.name in names], dtype=np.int32)
        idx_map.append(idx)

    return idx_map

def _centers_from_xyz(xyz_angstrom, idx_map):
    n = len(idx_map)
    centers = np.full((n, 3), np.nan, dtype=np.float32)
    for i, idxs in enumerate(idx_map):
        if idxs.size:
            centers[i] = xyz_angstrom[idxs].mean(axis=0, dtype=np.float32)
    return centers

def get_residue_distance_for_frame(trj: md.Trajectory,
                                    frame: int,
                                    atoms_spec = None,
                                    write_output: bool = True,
                                    idx_map: list[np.ndarray] | None = None) -> np.ndarray:
    """
    Docstring for get_residue_distance_for_frame
    
    :param trj: Description
    :type trj: md.Trajectory
    :param frame: Description
    :type frame: int
    :param atoms_spec: Description
    :param write_output: Description
    :type write_output: bool
    :param idx_map: Description
    :type idx_map: list[np.ndarray] | None
    :return: Description
    :rtype: ndarray[Any, Any]
    """
    
    tf = trj[frame - 1]
    if idx_map is None:
        idx_map = _build_ring_index_map(tf.topology, atoms_spec)

    if write_output:
        # very light progress only here
        sys.stdout.write(f"\rCenters for frame {frame}..."); sys.stdout.flush()

    xyzA = tf.xyz[0] * 10.0  # to Å
    centers = _centers_from_xyz(xyzA, idx_map)

    if write_output:
        sys.stdout.write(" distances..."); sys.stdout.flush()

    ssf = _pairwise_dist_numba(centers)

    if write_output:
        print(" done.")
    return ssf

# Globals each worker will fill once
_G = {"trj": None, "idx_map": None}

def _init_worker(trj_arg, atoms_spec):
    trj = md.load(trj_arg) if isinstance(trj_arg, str) else trj_arg
    _G["trj"] = trj
    _G["idx_map"] = _build_ring_index_map(trj.topology, atoms_spec)
    # keep any internal libraries single-threaded in each proc
    import os
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    try:
        import numba
        numba.set_num_threads(1)
    except Exception:
        pass

def _frame_worker(frame_1idx):
    tf = _G["trj"][frame_1idx - 1]
    xyzA = tf.xyz[0] * 10.0
    centers = _centers_from_xyz(xyzA, _G["idx_map"])
    return _pairwise_dist_numba(centers)

def get_residue_distance_for_trajectory(trj: md.Trajectory, 
                                         frames: typing.ArrayLike | str | set = {},
                                         atoms_spec = None,
                                         threads: int = 1,
                                         write_output: bool = True) -> np.ndarray:
    frames = SmartIndexingAction.parse_smart_index(frames) or list(range(1, trj.n_frames + 1))

    if threads == 1:
        idx_map = _build_ring_index_map(trj.topology, atoms_spec)
        out = []
        for k, f in enumerate(frames, 1):
            out.append(get_residue_distance_for_frame(trj, frame=f, atoms_spec=atoms_spec,
                                                       write_output=False, idx_map=idx_map))
            if write_output and (k % 25 == 0 or k == len(frames)):
                print(f"Frames done: {k}/{len(frames)}")
        return np.stack(out)

    if write_output:
        print("Note: suppressing per-frame prints in multiprocessing.")

    # Strongly prefer passing a path here if the traj is large.
    trj_arg = trj  

    # use SPAWN to avoid inheriting OpenMP state
    ctx = mp.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=threads,
        mp_context=ctx,
        initializer=_init_worker,
        initargs=(trj_arg, atoms_spec),
    ) as ex:
        out = list(ex.map(_frame_worker, frames))
    return np.stack(out)

@functools.wraps(get_residue_distance_for_trajectory)
def system_stacking_fingerprints(*args, **kwargs):
    return get_residue_distance_for_trajectory(*args, **kwargs) 

system_stacking_fingerprints.__doc__ = f"""
Alias for `get_residue_distance_for_trajectory()`.

{get_residue_distance_for_trajectory.__doc__}
"""

def get_frame_average(ssfs : typing.ArrayLike) -> typing.ArrayLike:
    '''
    Calculates an average System Stacking Fingerprint (SSF) across multiple SSFs

    Used to calculate an average SSF across multiple frames of a trajectory. Can
    average the result of `get_residue_distance_for_trajectory`

    Parameters
    ----------
    ssfs : numpy.typing.ArrayLike
        List or array of 2D NumPy arrays representing a pairwise distance matrix
        of an MD structure. All 2D NumPy arrays must be of the same dimenstions.
        Output of ``get_residue_distance_for_trajectory()``
        
    Returns
    -------
    avg_frame : numpy.typing.ArrayLike
        A single 2D NumPy array representing a pairwise distance matrix where each
        position i,j is the average distance from residue i to j across all matrices
        in frames.

    See Also
    --------
    get_residue_distance_for_trajectory : Calculates System Stacking Fingerprints (SSFs) for all residues across all frames of a trajectory
    system_stacking_fingerprints : Alias for get_residue_distance_for_trajectory

    Examples
    --------
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             residues = '2-5,13-16,23-31,46-51,65-76,88-104,122-141,164-175,184-198,288-289,401-415,420-430', 
    ...                             atoms = {'C2','C4','C6'})
    >>> ssfs = st.get_residue_distance_for_trajectory(filtered_traj, frames = '1-3', write_output = False)
    >>> avg_ssf = st.get_frame_average(ssfs)
    >>> avg_ssf.shape
    [FILL]
    '''
    avg_frame = np.mean(ssfs, axis = 0)
    return avg_frame 

def classify_residue_name(resname: str) -> str:
    """
    Classify a residue name into a coarse class:
      - 'nt'      : nucleotide
      - 'aa_arom' : aromatic amino acid (PHE, TYR, TRP, HIS variants)
      - 'aa_arg'  : arginine
      - 'other'   : anything else
    """
    if resname in _NUCLEOTIDE_NAMES:
        return "nt"
    if resname in _AROMATIC_AA_NAMES:
        return "aa_arom"
    if resname in _ARG_AA_NAMES:
        return "aa_arg"
    return "other"

def pair_type_from_classes(cls1: str, cls2: str) -> str:
    """
    Given two residue classes (from classify_residue_name), return a pair type:
      - 'nt-nt'
      - 'aa-aa'
      - 'aa-nt' (or 'nt-aa', symmetrized)
    """
    nt1 = (cls1 == "nt")
    nt2 = (cls2 == "nt")
    aa1 = cls1.startswith("aa")
    aa2 = cls2.startswith("aa")

    if nt1 and nt2:
        return "nt-nt"
    if aa1 and aa2:
        return "aa-aa"
    return "aa-nt"

def get_top_stacking(trj : md.Trajectory, matrix : typing.ArrayLike, csv : str = '',
                     n_events : int = 5, include_adjacent : bool = False) -> None:
    '''
    Returns top stacking residue pairs for a given System Stacking Fingerprint (SSF)

    Given a trajectory and a SSF made from `get_residue_distance_for_frame()` or `get_frame_average()`
    prints the residue pairings with the strongest stacking events (ie. the residue pairings
    with center of geometry distance closest to 3.5Å). 

    Parameters
    ----------    
    trj : md.Trajectory
        trajectory used to get the stacking fingerprint
    matrix : typing.ArrayLike
        Single-frame SSF created by ``get_residue_distance_for_frame()`` or ``get_frame_average()``
    csv : str, default = '',
        output filename of the tab-separated txt file to write data to. If empty, data printed to standard output
    n_events : int, default = 5
        maximum number of stacking events to display, if -1 display all residue pairings
    include_adjacent : bool, default = False
        True if adjacent residues should be included in the printed output

    See Also
    --------
    get_residue_distance_for_frame : Calculates System Stacking Fingerprint (SSF) between all residues in a given frame.
    get_residue_distance_for_trajectory : Calculates System Stacking Fingerprints (SSFs) for all residues across all frames of a trajectory
    system_stacking_fingerprints : Alias for get_residue_distance_for_trajectory

    Examples
    --------
    We can calculate the stacking events for a single frame:

    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             atoms = {'C2','C4','C6'})
    >>> ssf = st.get_residue_distance_for_frame(filtered_traj, frame = 2, write_output = False)
    >>> ssf.shape
    (252,252)
    >>> st.get_top_stacking(filtered_traj, ssf)
    Row     Column  Value
    197     195     3.50
    420     413     3.51
    94      127     3.51
    93      130     3.53
    117     108     3.38

    Or we can get most residue pairs that had the most stacking across many frames of a trajectory:

    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             atoms = {'C2','C4','C6'})
    >>> ssfs = st.get_residue_distance_for_trajectory(filtered_traj, frames = '1-3', write_output = False)
    >>> avg_ssf = st.get_frame_average(ssfs)
    >>> avg_ssf.shape
    (252, 252)
    >>> st.get_top_stacking(filtered_traj, ssf)
    Row     Column  Value
    130     93      3.56
    108     117     3.44
    195     197     3.61
    127     94      3.65
    47      167     3.66

    '''
    top_stacking_indices = np.argsort(np.abs(matrix - 3.5), axis = None)
    rows, cols = np.unravel_index(top_stacking_indices, matrix.shape)
    closest_values = matrix[rows, cols]

    if include_adjacent:
        # non_adjacent_indices includes adjacent indices in this case
        non_adjacent_indices = [(row, col, value) for row, col, value in zip(rows, cols, closest_values) if abs(row - col) > 0]
    else:
        non_adjacent_indices = [(row, col, value) for row, col, value in zip(rows, cols, closest_values) if abs(row - col) > 1]

    no_mirrored_indices = [] # keep only one side of x=y line, since mat[i,j] = mat[j,i]
    for row, col, value in non_adjacent_indices:
        if (col, row, value) not in no_mirrored_indices:
            no_mirrored_indices += [(row, col, value)]
    if n_events == -1: n_events = len(no_mirrored_indices) 
    no_mirrored_indices = no_mirrored_indices[:n_events]

    if csv:
        with open(csv, 'w') as csv_file:
            # New richer header
            csv_file.write(
                'Res1\tRes1_name\tRes1_class\t'
                'Res2\tRes2_name\tRes2_class\t'
                'Pair_type\tAvg_Dist\n'
            )
            for row, col, value in no_mirrored_indices:
                res1_res = trj.topology.residue(row)
                res2_res = trj.topology.residue(col)

                # 1-indexed residue IDs (PDB numbering)
                res1 = increment_residue(str(res1_res.resSeq))
                res2 = increment_residue(str(res2_res.resSeq))

                res1_name = res1_res.name
                res2_name = res2_res.name

                res1_class = classify_residue_name(res1_name)
                res2_class = classify_residue_name(res2_name)

                pair_type = pair_type_from_classes(res1_class, res2_class)

                csv_file.write(
                    f"{res1}\t{res1_name}\t{res1_class}\t"
                    f"{res2}\t{res2_name}\t{res2_class}\t"
                    f"{pair_type}\t{value:.2f}\n"
                )
    else:
        print(
            'Res1\tRes1_name\tRes1_class\t'
            'Res2\tRes2_name\tRes2_class\t'
            'Pair_type\tAvg_Dist'
        )
        for row, col, value in no_mirrored_indices:
            res1_res = trj.topology.residue(row)
            res2_res = trj.topology.residue(col)

            res1 = increment_residue(str(res1_res.resSeq))
            res2 = increment_residue(str(res2_res.resSeq))

            res1_name = res1_res.name
            res2_name = res2_res.name

            res1_class = classify_residue_name(res1_name)
            res2_class = classify_residue_name(res2_name)

            pair_type = pair_type_from_classes(res1_class, res2_class)

            print(
                f"{res1}\t{res1_name}\t{res1_class}\t"
                f"{res2}\t{res2_name}\t{res2_class}\t"
                f"{pair_type}\t{value:.2f}"
            )
 
def increment_residue(residue : str) -> str:
    '''
    Increments residue ID by 1
    
    Useful when converting from mdtraj 0-index residue naming to 1-indexed
    
    Parameters
    ----------
    residue : str
        The residue id given by trajectory.topology.residue(i)

    Returns
    -------
    incremented_id : str
        The residue id with the sequence number increased by 1

    Examples
    --------
    >>> increment_residue('G43')
    'G44'

    '''
    letter_part = ''.join(filter(str.isalpha, residue))
    number_part = ''.join(filter(str.isdigit, residue))
    incremented_number = str(int(number_part) + 1)
    return letter_part + incremented_number

def load_ssfs(file : str) -> typing.ArrayLike:
    """
    Loads a list of SSFs created by ``stacker -s ssf -d OUTFILE``

    Loads a list of SSFs where each element is an SSF from a given
    frame. This is ``OUTFILE`` when running ``stacker -s ssf -d OUTFILE``
    and quickly provides saved SSF data rather than recalculating SSFs
    for a given trajectory.

    Parameters
    ----------
    file : str 
        outfile from ``stacker -s ssf -d OUTFILE``. 
        Must be ``.txt`` or ``.txt.gz``.

    Returns
    -------
    ssfs : numpy.typing.ArrayLike
        List or array of 2D NumPy arrays representing a pairwise distance matrix
        of an MD structure. All 2D NumPy arrays must be of the same dimenstions.
        Output of ``get_residue_distance_for_trajectory()``

    See Also
    --------
    system_stacking_fingerprints : calculates ``ssfs`` rather than loading it like this function.
    get_residue_distance_for_trajectory : Alias for :func:`system_stacking_fingerprints()`
    get_frame_average : calculates average SSF for a trajectory. ``matrix`` parameter is the output of this
    """
    flattened_ssf = np.loadtxt(file)
    ssfs = flattened_ssf.reshape(flattened_ssf.shape[0], math.isqrt(flattened_ssf.shape[1]), math.isqrt(flattened_ssf.shape[1]))
    return ssfs

class MultiFrameTraj(Exception):
    """
    A multi-frame trajectory is passed to a one-frame function

    Raised if a multi-frame trajectory is passed to a function that
    only works on one trajectory (eg. calculate_residue_distance_vector())
    """
    pass

if __name__ == "__main__":
    trajectory_file = '../testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    topology_file = '../testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
    # Load test trajectory and topology
    trj = md.load(trajectory_file, top = topology_file)

    # "Correct" residue distances determined using PyMOL, a standard interface
    # for visualizing 3D molecules (distances limited to 3 decimal places)

    # calculate_residue_distance() tests
    tolerance = 1e-6
    assert round(calculate_residue_distance(trj[0], 426, 427).magnitude(), 3) - 7.525 < tolerance
    assert (round(calculate_residue_distance(trj[0], 3, 430).magnitude(), 3) - 22.043 < tolerance)
    ### Multi-frame exception
    try:
        round(calculate_residue_distance(trj[0:10], 3, 430).magnitude(), 3) - 22.043 < tolerance
    except MultiFrameTraj:
        print("MultiFrameTraj: calculate_residue_distance_vector() fails on multiple-frame trajectory")

    # create_axis_labels() test
    assert(create_axis_labels([0,1,2,3,4,5,6,7,8,9,10,11,12,98,99,100]) == ([0,10,12,13,15], [0,10,12,98,100]))
    assert(create_axis_labels([94,95,96,97,98,99,100,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428]) == ([0,6,7,17,27], [94,100,408,418,428]))
    ### No passed in residues exception
    try:
        assert(create_axis_labels([]) == ([],[]))
    except NoResidues:
        print("NoResidues: create_axis_labels() fails on empty residue list")

    # get_residue_distance_for_frame() test
    trj_three_residues = trj.atom_slice(trj.top.select('resi 407 or resi 425 or resi 426'))
    assert(np.all(np.vectorize(round)(get_residue_distance_for_frame(trj_three_residues, 2), 3) == np.array([[0,      8.231,   11.712], 
                                                                                                              [8.231,  0,       6.885], 
                                                                                                               [11.712, 6.885,   0]])))

    # display_arrays_as_video() tests
    residue_selection_query = 'resi 90 to 215'
    frames_to_include = [1,2,3,4,5]

    trj_sub = trj.atom_slice(trj.top.select(residue_selection_query))
    resSeqs = [res.resSeq for res in trj_sub.topology.residues]
    frames = get_residue_distance_for_trajectory(trj_sub, frames_to_include, threads = 5)
    get_top_stacking(trj_sub, frames[0])
    display_arrays_as_video([get_frame_average(frames)], resSeqs, seconds_per_frame=10)

    display_arrays_as_video(frames, resSeqs, seconds_per_frame=10)

    # All Residues one large matrix
    resSeqs = [res.resSeq for res in trj.topology.residues]
    print('\n')
    frames = [get_residue_distance_for_frame(trj, i) for i in range(1,2)]
    display_arrays_as_video(frames, resSeqs, seconds_per_frame=10, tick_distance=20)