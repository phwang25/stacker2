"""
Create and Analyze Pairwise Stacking Fingerprints (PSFs)

This module contains the functions to create and analyze Pairwise
Stacking Fingerprints. It allows the user to generate Polar Scatterplots
and Heatmaps of one residue movement relative to the other.

This pipeline takes a pdb file with multiple frames and two residues. 
It then calculates the r, rho, and theta values at each frame
as defined in the Bottaro paper (https://doi.org/10.1093/nar/gku972), 
and exports these values to a .csv file.
"""

from __future__ import print_function
import math
import csv
import mdtraj as md
from .vector import *
from .file_manipulation import filter_traj_to_pdb
from .visualization import create_parent_directories
import os
import functools
from typing import Optional, Set, Tuple

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

def collect_atom_locations_by_frame(trj: md.Trajectory, residue: int, atom: str) -> list:
    """
    Creates a list of all atom locations for a particular atom and residue number per frame.

    Curates a list coords_by_frame where coords_by_frame[i] is the (x, y, z) positions 
    of a provided ``atom`` in a residue ``residue`` at the ith frame.

    Parameters
    ----------
    trj : md.Trajectory
        Trajectory to analyze.
    residue : int
        The 0-indexed residue number of the residue where `atom_id` is found (PDB Column 5).
    atom : str
        The name of the atom to get coordinates for (PDB Column 2).

    Returns
    -------
    coords_by_frame : list
        List of (x, y, z) coordinates of `atom_id` in `residue_num` for each frame.

    Notes
    -----
    `residue_num` must be 0-indexed to match how mdtraj.Trajectory indexes residues.

    See Also
    --------
    Base : Python Class that represents a nucleotide base
    calculate_bottaro_values_for_frame : Calculates the r, rho, and theta values as expressed in the Bottaro paper.
    create_base_from_coords_list : Combines C2, C4, C6 positions with midpoint positions for a given frame
    
    Examples
    --------
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                              'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                              atoms = {'C2','C4','C6'})
    >>> >>> st.collect_atom_locations_by_frame(filtered_traj, residue = 3, atom = "C2")
    [(58.794, 59.636, 49.695), (59.185005, 58.797, 50.137), (59.379, 58.553005, 49.853), 
    (58.76, 59.068, 49.681), (59.003, 59.054, 49.878002), (59.049, 58.967, 50.051), 
    (59.219, 58.476006, 49.948), (58.948, 58.588005, 50.085), (58.922, 58.747, 49.766003), 
    (59.124, 58.916, 49.978004)]

    """
    topology = trj.topology
    number_of_frames = trj.n_frames
    atomic_index = topology.select("name " + atom + " and residue " + str(residue))[0]

    # multiply by 10 to convert nanometer units in trajectory.xyz to Angstroms
    coords_by_frame = [tuple(trj.xyz[frame_idx, atomic_index,:] * 10) for frame_idx in range(0, number_of_frames)]
    return coords_by_frame

def calc_center_3pts(a: Vector, b: Vector, c: Vector) -> Vector:
    """
    Finds the average x, y, z position of three (x, y, z) Vectors.

    Takes in three Vectors generated using the `pdb.xyz` method and finds their center. Works 
    with three points that make up a triangle (like those within a 6-member ring).

    Parameters
    ----------
    a : Vector
        x, y, z coordinates of the first point.
    b : Vector
        x, y, z coordinates of the second point.
    c : Vector
        x, y, z coordinates of the third point.

    Returns
    -------
    midpoint : Vector
        One Vector with (x, y, z) coordinates at the center of the three input Vectors.

    See Also
    --------
    Base : Python Class that represents a nucleotide base
    calculate_bottaro_values_for_frame : Calculates the r, rho, and theta values as expressed in the Bottaro paper.
    create_base_from_coords_list : Combines C2, C4, C6 positions with midpoint positions for a given frame

    Examples
    --------
    >>> import stacker as st
    >>> print(st.calc_center_3pts(st.Vector(0,0,0), st.Vector(1,1,1), st.Vector(2,2,2)))
    [ 1.0
      1.0
      1.0 ]

    """
    vectorized = (a.components + b.components + c.components ) / 3
    midpoint = Vector(*vectorized)
    return midpoint

class Base:
    """
    Represents a nucleotide base with x, y, z coordinates for C2, C4, and C6 atoms,
    and their average position at a single frame.

    This class defines a data type 'Base' that consists of coordinates for the atoms
    C2, C4, and C6, as well as the midpoint of these atoms for a single residue at a single frame.

    Attributes
    ----------
    c2_coords : Vector
        (x, y, z) coordinates for the atom C2.
    c4_coords : Vector
        (x, y, z) coordinates for the atom C4.
    c6_coords : Vector
        (x, y, z) coordinates for the atom C6.
    midpoint_coords : Vector
        (x, y, z) coordinates representing the midpoint of C2, C4, and C6.

    """
    def __init__(self, c2_coords: tuple, c4_coords: tuple, c6_coords: tuple, midpoint_coords: tuple) -> None:
        """
        Initialize a Base instance.

        Parameters
        ----------
        c2_coords : tuple
            (x, y, z) coordinates for C2.
        c4_coords : tuple
            (x, y, z) coordinates for C4.
        c6_coords : tuple
            (x, y, z) coordinates for C6.
        midpoint_coords : tuple
            (x, y, z) coordinates for the midpoint.

        """
        self.c2_coords = Vector(*c2_coords)
        self.c4_coords = Vector(*c4_coords)
        self.c6_coords = Vector(*c6_coords)
        self.midpoint_coords = Vector(*midpoint_coords)

def create_base_from_coords_list(frame: int, C2_coords: list, C4_coords: list, C6_coords: list, midpoint_coords: list) -> Base:
    """
    Combines C2, C4, C6 positions with midpoint positions for a given frame.

    Takes a frame (0-indexed) and outputs a Base instance with the x, y, z locations of C2, C4, C6, and midpoint 
    of the same residue at that frame. 

    Parameters
    ----------
    frame : int
        Frame number (0-indexed).
    C2_coords : list
        List of (x, y, z) coordinates of C2 atom in some residue for each frame.
    C4_coords : list
        List of (x, y, z) coordinates of C4 atom in some residue for each frame.
    C6_coords : list
        List of (x, y, z) coordinates of C6 atom in some residue for each frame.
    midpoint_coords : list
        List of (x, y, z) coordinates of the midpoint of residue for each frame.

    Returns
    -------
    Base
        An instance of Base with coordinates at the specified frame.

    See Also
    --------
    calc_center_3pts : Finds the average x, y, z position of three (x, y, z) Vectors

    Notes
    -----
    `frame` is 0-indexed because this is how `mdtraj` handles mdtraj.Trajectory
    
    """
    midpoint_tuple = (midpoint_coords[frame].x, midpoint_coords[frame].y, midpoint_coords[frame].z)
    return Base(C2_coords[frame], C4_coords[frame], C6_coords[frame], midpoint_tuple)

def correct_theta_sign(rho: Vector, y_axis: Vector, theta: float) -> float:
    """
    Corrects the sign of an angle theta with the x-axis within a plane defined by a given y-axis.

    When calculating the angle theta between two vectors in 3D space, once the vectors move
    >180 degrees apart, the angle becomes the shortest path. To have 360 degrees of freedom, we calculate
    theta within the plane by checking if it faces the same direction as a vector y, and correcting
    otherwise.

    Parameters
    ----------
    rho : Vector
        The vector compared to the x-axis to form theta.
    y_axis : Vector
        Directional vector to define a plane with x-axis; orthogonal to x-axis.
    theta : float
        The calculated angle of rho with x-axis to be corrected.

    Returns
    -------
    float
        Theta as calculated on the plane.

    """
    proj_rho_on_y = rho.calculate_projection(y_axis)
    if y_axis.x != 0:
        opposite_direction = (proj_rho_on_y.x / y_axis.x < 0)
        if opposite_direction:
            theta = 360 - theta
    return theta

def calculate_bottaro_values_for_frame(perspective_base_coords: Base, viewed_midpoint: Vector) -> list:
    """
    Calculates the r, rho, and theta values as expressed in the Bottaro paper.

    Calculates the r, rho, and theta values between two nucleotides in a single frame
    as presented in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972).

    Parameters
    ----------
    perspective_base_coords : Base
        List of the x, y, z coords of C2, C4, C6 and their midpoint for 
        the perspective residue in a single frame.
    viewed_midpoint : Vector
        x, y, z position of the midpoint of the viewed nucleotide at the same frame.

    Returns
    -------
    list
        A list containing 3 floats from Bottaro; structure: [r_dist, rho_dist, theta].

    References
    ----------
    [1] Sandro Bottaro, Francesco Di Palma, Giovanni Bussi, The role of nucleobase interactions in 
    RNA structure and dynamics, Nucleic Acids Research, Volume 42, Issue 21, 1 December 2014, 
    Pages 13306–13314, https://doi.org/10.1093/nar/gku972

    """
    r_vector = viewed_midpoint - perspective_base_coords.midpoint_coords
    r_magnitude = r_vector.magnitude()
    
    # 2 Vectors define the plane of the perspective nucleotide
    x_axis = perspective_base_coords.c2_coords - perspective_base_coords.midpoint_coords
    vector_midpoint_to_C4 = perspective_base_coords.c4_coords - perspective_base_coords.midpoint_coords

    normal_vector_to_plane = x_axis.calculate_cross_product(vector_midpoint_to_C4)
    
    # y-axis inside plane to get correct theta value in 3D
    y_axis = x_axis.calculate_cross_product(normal_vector_to_plane)

    proj_r_on_normal = r_vector.calculate_projection(normal_vector_to_plane)

    rho = r_vector - proj_r_on_normal
    rho_dist = rho.magnitude()
    
    x_axis_dot_rho = x_axis.x * rho.x + x_axis.y * rho.y + x_axis.z * rho.z
    denominator = x_axis.magnitude() * rho_dist
    cos_theta = x_axis_dot_rho / denominator

    # edge cases where rounding leads to minor error
    if cos_theta > 1: cos_theta = 1 
    if cos_theta < -1: cos_theta = -1

    theta = math.degrees(math.acos(cos_theta))
    corrected_theta = correct_theta_sign(rho, y_axis, theta)
    
    values = [r_magnitude, rho_dist, corrected_theta]
    return values

def write_bottaro_to_csv(pdb: str = '', 
                         outcsv: str = '', 
                         pers_res: int = -1, 
                         view_res: int = -1,
                         res1_atoms: set = {"C2", "C4", "C6"}, 
                         res2_atoms: set = {"C2", "C4", "C6"}, 
                         index: int = 1) -> None:
    """
    Write the r, rho, and theta values used to make a Pairwise Stacking Fingerprint (PSF)
    from a trajectory PDB to a CSV.

    Calculates the r, rho, and theta values as described in Bottaro et al. from a
    perspective nucleotide residue to a viewed nucleotide residue per frame. Writes the 
    results to a CSV file.

    Parameters
    ----------
    pdb : str
        Filename of PDB containing information for ONLY two residues (perspective and viewed
        nucleotide) at each frame.
    outcsv : str
        Filename of CSV file to write to.
    pers_res : int, default = -1
        Residue index of the perspective residue whose plane to project onto (0-/1-index changed by
        `index` variable, default 1-indexed). If -1, a 2-residue PDB is assumed and perspective id is
        the first res_id.
    view_res : int, default = -1
        Residue index of the viewed residue whose midpoint to project to perspective residue plane 
        (0-/1-index changed by index variable, default 1-indexed). If -1, a 2-residue PDB is assumed 
        and viewed id is the second res_id.
    res1_atoms : set, default = {"C2", "C4", "C6"}
        Set of the atom names (e.g., "C2", "C4", "C6") to use from
        residue 1 to find center of geometry for perspective nucleotide.
    res2_atoms : set, default = {"C2", "C4", "C6"}
        Set of the atom names (e.g., "C2", "C4", "C6") to use from
        residue 2 to find center of geometry for viewed nucleotide.
    index : int, default = 1
        Index of the residues. 1-indexed (default) means residue ids start at 1.
        cpptraj uses 1-indexed residues. mdtraj PDB outputs will be 0-indexed.
        
    See Also
    --------
    filter_traj_to_pdb : convert a trajectory and topology file into a single filtered PDB file to input here
    visualize_two_residue_movement_scatterplot : Visualize Output CSV data as a PSF scatterplot from this data
    visualize_two_residue_movement_heatmap : Visualize Output CSV data as a PSF heatmap from this data
        
    Examples
    --------
    >>> import stacker as st
    >>> trajectory_file = 'testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    >>> topology_file = 'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
    >>> pdb_filename = 'testing/script_tests/residue_movement/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb'
    >>> output_csv_name = "testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv"
    >>> perspective_residue = 426 # 1-indexed
    >>> viewed_residue = 427 # 1-indexed
    >>> st.filter_traj_to_pdb(trj_file=trajectory_file, top_file=topology_file, pdb=pdb_filename,
    ...                        residues={perspective_residue,viewed_residue}, atoms={"C2", "C4", "C6"})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    WARNING: Output file atom, residue, and chain indices are zero-indexed
    Filtered trajectory written to:  testing/script_tests/residue_movement/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb
    >>> st.write_bottaro_to_csv(pdb_filename, output_csv_name, pers_res=perspective_residue, view_res=viewed_residue)
    Output values written to testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv
    >>> print("".join(open(output_csv_name).readlines()[:10]))
    frame,r_dist,rho_dist,theta
    0,7.5253415,6.5321836,204.02934901525177
    1,6.884639,6.0513134,199.40647902703924
    2,7.301847,6.151191,205.1906453260924
    3,6.5815425,5.461494,199.5421877249345
    4,7.0760417,5.3919506,204.0150121540755
    5,7.2589145,6.3483577,201.32674968542617
    6,7.4929285,6.414151,205.92967194025135
    7,7.1484976,6.035165,202.88441276229827
    8,7.344863,5.541237,217.30043061558888

    References
    ----------
    [1] Sandro Bottaro, Francesco Di Palma, Giovanni Bussi, The role of nucleobase interactions in 
    RNA structure and dynamics, Nucleic Acids Research, Volume 42, Issue 21, 1 December 2014, 
    Pages 13306–13314, https://doi.org/10.1093/nar/gku972

    """
    # Keep atomname order consistent between runs
    res1_atoms = list(res1_atoms)
    res2_atoms = list(res2_atoms)
    res1_atoms.sort()
    res2_atoms.sort()

    res1_atom1,res1_atom2,res1_atom3 = res1_atoms
    res2_atom1,res2_atom2,res2_atom3 = res2_atoms

    pdb = md.load(pdb)
    number_of_frames = pdb.n_frames

    if index == 1: # correct for 0-index res_id of mdtraj
        pers_res -= 1
        view_res -= 1

    if pers_res == -1 or pers_res == -2:
        topology = pdb.topology
        pers_res = [residue for residue in topology.residues][0].resSeq

    if view_res == -1 or view_res == -2:
        topology = pdb.topology
        view_res = [residue for residue in topology.residues][1].resSeq

    residue1_C2_list = collect_atom_locations_by_frame(pdb, pers_res, res1_atom1)
    residue1_C4_list = collect_atom_locations_by_frame(pdb, pers_res, res1_atom2)
    residue1_C6_list = collect_atom_locations_by_frame(pdb, pers_res, res1_atom3)
    residue2_C2_list = collect_atom_locations_by_frame(pdb, view_res, res2_atom1)
    residue2_C4_list = collect_atom_locations_by_frame(pdb, view_res, res2_atom2)
    residue2_C6_list = collect_atom_locations_by_frame(pdb, view_res, res2_atom3)
    
    residue1_midpoint_list = [calc_center_3pts(Vector(*residue1_C2_list[i]),
                                               Vector(*residue1_C4_list[i]),
                                               Vector(*residue1_C6_list[i])) for i in range(0, number_of_frames)]
    
    residue2_midpoint_list = [calc_center_3pts(Vector(*residue2_C2_list[i]),
                                               Vector(*residue2_C4_list[i]),
                                               Vector(*residue2_C6_list[i])) for i in range(0, number_of_frames)]
    
    fields = ['frame','r_dist','rho_dist', 'theta']
    rows=[]
    for i in range(0,number_of_frames):
        residue1_base = create_base_from_coords_list(i, residue1_C2_list, residue1_C4_list, residue1_C6_list, residue1_midpoint_list)
        frame_values = calculate_bottaro_values_for_frame(residue1_base,residue2_midpoint_list[i])
        row = [i]+frame_values
        rows.append(row)
    
    filename = outcsv
    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile) 
        csvwriter.writerow(fields) 
        csvwriter.writerows(rows)
    print("Output values written to " + outcsv)

# ————————————————————————————————————————————————————————————————————————————————
# NEW VERSION BELOW

def determine_default_ring_atoms(res_name: str) -> set:
    """Return default atom names defining the planar aromatic group for PSF v2."""
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

def determine_default_anchors(res_name: str) -> Tuple[str, str]:
    """Return (A0, A1) anchor atom names for SVD-based PSF orientation.

    A0 defines the theta=0 direction in the perspective plane.
    A1 fixes handedness/sign of the SVD normal (prevents theta mirroring across frames).
    """
    # Nucleotides: preserve legacy Bottaro convention
    if res_name in _NUCLEOTIDE_NAMES:
        return ("C2", "C4")

    # Aromatic AAs: attachment point + deterministic neighbor
    if res_name in _PHE_NAMES or res_name in _TYR_NAMES or res_name in _TRP_NAMES:
        return ("CG", "CD1")
    if res_name in _HIS_NAMES:
        return ("CG", "ND1")

    # Guanidinium: attachment-side NE + NH1 to fix handedness
    if res_name in _ARG_NAMES:
        return ("NE", "NH1")

    raise ValueError(f"No default anchors defined for residue name '{res_name}'.")

def _select_atom_index(topology: md.Topology, residue_idx: int, atom_name: str) -> int:
    sel = topology.select(f"name {atom_name} and residue {residue_idx}")
    if sel.size == 0:
        raise ValueError(f"Could not find atom '{atom_name}' in residue index {residue_idx}.")
    return int(sel[0])

def _select_atom_indices(topology: md.Topology, residue_idx: int, atom_names: Set[str]) -> np.ndarray:
    idx = []
    for a in topology.residue(residue_idx).atoms:
        if a.name in atom_names:
            idx.append(a.index)
    return np.array(idx, dtype=np.int32)

def _svd_plane_from_coords(coordsA: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (center, unit normal) for best-fit plane from coordsA (Å)."""
    C = coordsA.mean(axis=0)
    Y = coordsA - C
    # SVD on N x 3; normal is smallest right-singular vector
    _, _, Vt = np.linalg.svd(Y, full_matrices=False)
    n = Vt[-1]
    n_norm = float(np.linalg.norm(n))
    if n_norm == 0.0:
        raise ValueError("SVD produced a zero normal vector (degenerate coordinates).")
    n = n / n_norm
    return C, n

def _fix_normal_sign(n_svd: np.ndarray, C: np.ndarray, A0: np.ndarray, A1: np.ndarray) -> np.ndarray:
    """Flip n_svd if needed so it matches anchor-handedness from (A0, A1)."""
    v1 = A0 - C
    v2 = A1 - C
    cross = np.cross(v1, v2)
    cross_norm = float(np.linalg.norm(cross))
    if cross_norm == 0.0:
        # anchors colinear w/ centroid; can't define handedness -> do nothing
        return n_svd
    n_anchor = cross / cross_norm
    if float(np.dot(n_svd, n_anchor)) < 0.0:
        return -n_svd
    return n_svd

def _in_plane_basis(C: np.ndarray, n: np.ndarray, A0: np.ndarray, ring_coordsA: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (xhat, yhat) in the plane using A0 as theta=0 reference."""
    u = A0 - C
    u_par = u - np.dot(u, n) * n
    u_norm = float(np.linalg.norm(u_par))

    if u_norm < 1e-8:
        # Fallback: choose the ring atom farthest from C after projection
        best = None
        best_norm = -1.0
        for p in ring_coordsA:
            v = p - C
            v_par = v - np.dot(v, n) * n
            vn = float(np.linalg.norm(v_par))
            if vn > best_norm:
                best_norm = vn
                best = v_par
        if best is None or best_norm < 1e-8:
            raise ValueError("Could not construct in-plane x-axis (ring collapsed or all points on normal).")
        xhat = best / best_norm
    else:
        xhat = u_par / u_norm

    # yhat = np.cross(n, xhat) # right handedness (new method)
    yhat = np.cross(xhat, n)  # legacy Bottaro handedness (left handed)
    y_norm = float(np.linalg.norm(yhat))
    if y_norm == 0.0:
        raise ValueError("Could not construct in-plane y-axis (degenerate basis).")
    yhat = yhat / y_norm
    return xhat, yhat

def _bottaro_values_svd(pers_ringA: np.ndarray,
                        view_ringA: np.ndarray,
                        pers_A0A: np.ndarray,
                        pers_A1A: np.ndarray,
                        view_A0A: Optional[np.ndarray] = None,
                        view_A1A: Optional[np.ndarray] = None, ) -> Tuple[float, float, float, float]:
    
    """Compute (r_dist, rho_dist, theta_deg, tilt_deg) using SVD plane(s)."""
    # Perspective plane
    Ci, ni = _svd_plane_from_coords(pers_ringA)
    ni = _fix_normal_sign(ni, Ci, pers_A0A, pers_A1A)
    xhat, yhat = _in_plane_basis(Ci, ni, pers_A0A, pers_ringA)

    # Viewed center
    Cj = view_ringA.mean(axis=0)

    r = Cj - Ci
    r_dist = float(np.linalg.norm(r))

    rho_vec = r - np.dot(r, ni) * ni
    rho_dist = float(np.linalg.norm(rho_vec))

    # theta in [0, 360)
    x = float(np.dot(rho_vec, xhat))
    y = float(np.dot(rho_vec, yhat))
    theta = (math.degrees(math.atan2(y, x)) % 360.0)

    # Tilt angle between planes (optional but useful)
    tilt_deg = float('nan')
    if view_A0A is not None and view_A1A is not None:
        Cj2, nj = _svd_plane_from_coords(view_ringA)
        nj = _fix_normal_sign(nj, Cj2, view_A0A, view_A1A)
        # parallelness in [0, 90]
        dot = abs(float(np.dot(ni, nj)))
        dot = min(1.0, max(-1.0, dot))
        tilt_deg = math.degrees(math.acos(dot))

    return r_dist, rho_dist, theta, tilt_deg

def write_bottaro_to_csv_svd(pdb: str,
                             outcsv: str,
                             pers_res: int = -1,
                             view_res: int = -1, 
                             pers_ring_atoms: Optional[Set[str]] = None,
                             view_ring_atoms: Optional[Set[str]] = None,
                             index: int =1) -> None:
    """
    PSF v2 writer: compute Bottaro (r, rho, theta) using **full-ring COG + SVD plane**.

    Works for nucleotides (pur/pyr/INO), aromatic amino acids (PHE/TYR/TRP/HIS variants),
    and ARG guanidinium.

    Notes
    -----
    - `pers_res` / `view_res` are residue ids in the *input pdb*; if `index=1` they are treated as 1-indexed.
    - For easiest setup, generate the PDB with:
        `filter_traj_to_pdb(..., residues={pers_res, view_res}, atoms=set())`
      so all atoms for those residues are present.

    Output columns
    --------------
    frame,r_dist,rho_dist,theta,tilt_deg
    """
    trj = md.load(pdb)
    top = trj.topology

    if index == 1: # correct for 0-index res_id of mdtraj
        pers_res -= 1
        view_res -= 1
    
    if pers_res == -1 or pers_res == -2:
        pers_res = [residue for residue in top.residues][0].resSeq

    if view_res == -1 or view_res == -2:
        view_res = [residue for residue in top.residues][1].resSeq

    # Residue indices in topology (assuming PDB only has 2 residues)
    pers_index = 0
    view_index = 1

    if top.residue(pers_index).resSeq != pers_res or top.residue(view_index).resSeq != view_res:
        raise ValueError(
            f"PDB residue order mismatch: residue[0].resSeq={top.residue(0).resSeq}, "
            f"residue[1].resSeq={top.residue(1).resSeq}, but expected pers_res={pers_res}, view_res={view_res}."
        )

    pers_name = top.residue(pers_index).name
    view_name = top.residue(view_index).name

    if pers_ring_atoms is None:
        pers_ring_atoms = determine_default_ring_atoms(pers_name)
    if view_ring_atoms is None:
        view_ring_atoms = determine_default_ring_atoms(view_name)

    if not pers_ring_atoms:
        raise ValueError(f"No ring atom set found for perspective residue '{pers_name}'.")
    if not view_ring_atoms:
        raise ValueError(f"No ring atom set found for viewed residue '{view_name}'.")

    pers_A0_name, pers_A1_name = determine_default_anchors(pers_name)
    view_A0_name, view_A1_name = determine_default_anchors(view_name)

    pers_ring_idx = _select_atom_indices(top, pers_index, pers_ring_atoms)
    view_ring_idx = _select_atom_indices(top, view_index, view_ring_atoms)

    if pers_ring_idx.size < 3:
        raise ValueError(f"Perspective residue {pers_res+1} has <3 ring atoms in the PDB (found {pers_ring_idx.size}).")
    if view_ring_idx.size < 3:
        raise ValueError(f"Viewed residue {view_res+1} has <3 ring atoms in the PDB (found {view_ring_idx.size}).")

    pers_A0_idx = _select_atom_index(top, pers_res, pers_A0_name)
    pers_A1_idx = _select_atom_index(top, pers_res, pers_A1_name)
    view_A0_idx = _select_atom_index(top, view_res, view_A0_name)
    view_A1_idx = _select_atom_index(top, view_res, view_A1_name)

    fields = ['frame', 'r_dist', 'rho_dist', 'theta', 'tilt_deg']
    rows = []

    for f in range(trj.n_frames):
        xyzA = trj.xyz[f] * 10.0  # nm -> Å

        pers_ringA = xyzA[pers_ring_idx]
        view_ringA = xyzA[view_ring_idx]

        pers_A0A = xyzA[pers_A0_idx]
        pers_A1A = xyzA[pers_A1_idx]
        view_A0A = xyzA[view_A0_idx]
        view_A1A = xyzA[view_A1_idx]

        r_dist, rho_dist, theta, tilt_deg = _bottaro_values_svd(
            pers_ringA, view_ringA, pers_A0A, pers_A1A, view_A0A, view_A1A
        )

        rows.append([f, r_dist, rho_dist, theta, tilt_deg])

    with open(outcsv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)

    print(f"Output values written to {outcsv}")

@functools.wraps(write_bottaro_to_csv_svd)
def write_psf_data_svd(*args, **kwargs):
    """Alias for `write_bottaro_to_csv_svd()` (PSF v2: SVD plane + full-ring support)."""
    return write_bottaro_to_csv_svd(*args, **kwargs)

write_psf_data_svd.__doc__ = f"""
Alias for `write_bottaro_to_csv_svd()`.

{write_bottaro_to_csv_svd.__doc__}
"""

# ————————————————————————————————————————————————————————————————————————————————
# FOR OBTAINING ANGLE BETWEEN PLANES 
# #TODO (REMOVE LATER)

def ring_center_and_normal(coords: np.ndarray) -> tuple[Vector, Vector]:
    """
    Given an array of ring atom coordinates, fit a plane and return its center and normal.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n_atoms, 3) with x, y, z coordinates (any length units).

    Returns
    -------
    center : Vector
        Center of geometry of the ring.
    normal : Vector
        Unit normal vector to the best-fit plane through the atoms.
        If the fit is degenerate, returns Vector(0,0,0) as the normal.
    """
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coords must be of shape (n_atoms, 3)")

    # Center of geometry
    center_xyz = coords.mean(axis=0)
    centered = coords - center_xyz

    # Degenerate case: all points identical
    if np.allclose(centered, 0.0):
        return Vector(*center_xyz), Vector(0.0, 0.0, 0.0)

    # SVD: smallest singular vector is the normal
    # centered ~ U S V^T, rows are points, columns are x,y,z
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]  # vector corresponding to smallest singular value

    norm = np.linalg.norm(normal)
    if norm < 1e-8:
        normal = np.array([0.0, 0.0, 0.0])
    else:
        normal = normal / norm

    return Vector(*center_xyz), Vector(*normal)

def get_ring_atom_indices(topology: md.Topology, resseq: int, ring_atom_names: tuple[str, ...]) -> np.ndarray:
    """
    Return atom indices for the specified ring atoms in a given residue (by resSeq).

    Parameters
    ----------
    topology : md.Topology
    resseq : int
        PDB-style residue sequence number (resSeq).
    ring_atom_names : tuple of str
        Atom names that define the ring (e.g. ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ') for PHE).

    Returns
    -------
    indices : np.ndarray
        Array of atom indices.
    """
    name_clause = " or ".join(f"name {nm}" for nm in ring_atom_names)
    # Adjust resseq handling if your PDB is 0/1 off; this follows your existing resSeq usage
    selection_str = f"({name_clause}) and resSeq {resseq}"
    indices = topology.select(selection_str)
    return indices

def get_residue_ring_center_and_normal(trj: md.Trajectory, frame: int, resseq: int, ring_atom_names: tuple[str, ...], coords_in_angstrom: bool = True,) -> tuple[Vector, Vector]:
    """
    For a single residue in a single frame, compute center and normal of its aromatic ring.

    Parameters
    ----------
    trj : md.Trajectory
        Trajectory object.
    frame : int
        0-indexed frame index.
    resseq : int
        PDB-style resSeq number of the residue.
    ring_atom_names : tuple of str
        Atom names defining the ring.
    coords_in_angstrom : bool, optional
        If True, convert nm -> Å for the returned center.

    Returns
    -------
    center : Vector
        Center of the ring.
    normal : Vector
        Unit normal vector to the ring plane.
    """
    topology = trj.topology
    atom_indices = get_ring_atom_indices(topology, resseq, ring_atom_names)

    if atom_indices.size == 0:
        # No atoms found; return zeros
        return Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0)

    coords_nm = trj.xyz[frame, atom_indices, :]  # shape (n_atoms, 3), in nm
    if coords_in_angstrom:
        coords = coords_nm * 10.0  # nm -> Å
    else:
        coords = coords_nm

    return ring_center_and_normal(coords)

def angle_between_normals_deg(n1: Vector, n2: Vector) -> float:
    """
    Compute the angle in degrees between two normal vectors.

    Returns NaN if any normal is zero.
    """
    mag1 = n1.magnitude()
    mag2 = n2.magnitude()
    if mag1 == 0.0 or mag2 == 0.0:
        return float("nan")

    dot = n1.calculate_dot_product(n2) / (mag1 * mag2)
    # Clamp numerical noise
    dot = max(-1.0, min(1.0, dot))
    angle_rad = math.acos(dot)
    return math.degrees(angle_rad)

def distance_and_normal_angle_for_pair(trj: md.Trajectory, frame: int, resseq1: int, ring_atoms1: tuple[str, ...], resseq2: int, ring_atoms2: tuple[str, ...],) -> tuple[float, float]:
    """
    Compute COG distance and ring-normal angle (degrees) for a pair of residues in one frame.

    Returns
    -------
    distance : float
        Distance between ring centers (Å).
    angle_deg : float
        Angle between ring normals (0–180 degrees). NaN if degenerate.
    """
    c1, n1 = get_residue_ring_center_and_normal(trj, frame, resseq1, ring_atoms1)
    c2, n2 = get_residue_ring_center_and_normal(trj, frame, resseq2, ring_atoms2)

    # Distance between centers
    d_vec = c2 - c1
    distance = d_vec.magnitude()

    angle_deg = angle_between_normals_deg(n1, n2)
    return distance, angle_deg

# ————————————————————————————————————————————————————————————————————————————————

if __name__ == "__main__":
    trajectory_file = 'stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    topology_file = 'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
    output_csv_name = "stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv"
    perspective_residue = 426 # 1-indexed
    viewed_residue = 427 # 1-indexed
    create_parent_directories(output_csv_name)

    ########OPTIONAL VARS#######
    perspective_atom1_name = "C2"
    perspective_atom2_name = "C4"
    perspective_atom3_name = "C6"
    viewed_atom1_name = "C2"
    viewed_atom2_name = "C4"
    viewed_atom3_name = "C6"
    ############################

    pdb_filename = 'stacker/testing/script_tests/residue_movement/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb'
    filter_traj_to_pdb(trj_file=trajectory_file, top_file=topology_file, pdb=pdb_filename,
                       residues={perspective_residue,viewed_residue}, atoms={"C2", "C4", "C6"})

    # Two Residue movement test 10 frames
    write_bottaro_to_csv(pdb_filename, 
                         output_csv_name, pers_res=perspective_residue, view_res=viewed_residue,
                         res1_atoms={perspective_atom1_name, perspective_atom2_name, perspective_atom3_name}, 
                         res2_atoms={viewed_atom1_name,viewed_atom2_name,viewed_atom3_name})
    
    multiframe_pdb = 'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb'
    multiframe_csv = 'stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv'
    write_bottaro_to_csv(multiframe_pdb, multiframe_csv)