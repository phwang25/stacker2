"""Filter Trajectory Files

This module contains functions for manipulating trajectory
files. This includes filtering a trajectory to desired atoms,
converting trajectory filetype, and outputting Python trajectories
to other filetypes (eg. prmtop, mdcrd, pdb)
"""

import mdtraj as md
from numpy import typing
import argparse

def filter_traj(trj_file : str, top_file : str, 
                residues : str | set = {}, atoms : set = {}) -> md.Trajectory:
    '''
    Filters an input trajectory to only the specified atoms and residues

    Filteres an input trajectory that contains all of the atoms in a topology to only
    the desired atoms at the desired residues (eg. the atoms necessary to find the 
    center of geometry of a residue). If ``residues`` or ``atoms`` are
    empty, all residues or atoms are included respectively.

    Parameters
    ----------
    trj_file : str
        filepath of the trajectory
    top_file : str
        filepath of the topology of the molecule
    residues : set or str
        1-indexed residue numbers of residues to keep in the trajectory.
        Accepts smart-indexed str representing a list of residues (e.g '1-5,6,39-48').
        If Empty, include all residues.
    atoms : set 
        atomnames to keep in the trajectory. If Empty, include all atoms.
        
    Returns
    -------
    filtered_trajectory : mdtraj.Trajectory
        a trajectory object representing the filtered structure across all frames

    See Also
    --------
    filter_traj_to_pdb : Filters an input trajectory to only the specified 
                         atoms and residues and outputs to pdb
    mdtraj.Trajectory : The Trajectory object in mdtraj package
    
    Notes
    -----
    Inputed trajectory should have 1-indexed Residue Indices, 
    Outputed trajectory object will be 0-indexed.

    Examples
    --------
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             residues = {426,427}, 
    ...                             atoms = {'C2','C4','C6'})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    >>> table, bonds = filtered_traj.topology.to_dataframe()
    >>> print(table)
    serial name element  resSeq resName  chainID segmentID
    0   None   C6       C     425       G        0          
    1   None   C2       C     425       G        0          
    2   None   C4       C     425       G        0          
    3   None   C6       C     426       C        0          
    4   None   C4       C     426       C        0          
    5   None   C2       C     426       C        0       

    Residue Number support SmartIndexing
    
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('../testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             '../testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                              residues = '1-16,25,50-57', 
    ...                              atoms = {'C2','C4','C6'})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    >>> filtered_traj
    <mdtraj.Trajectory with 10 frames, 75 atoms, 25 residues, without unitcells at 0x1156c3ed0>

    '''
    print("WARNING: Residue Indices are expected to be 1-indexed")
    
    print("Reading trajectory...")
    trajectory = md.load(trj_file, top = top_file)
    
    print("Reading topology...")
    topology = trajectory.topology
    
    print("Filtering trajectory...")
    residues = SmartIndexingAction.parse_smart_index(residues)

    # make resSeq 0-indexed for mdtraj query
    residues = {resnum-1 for resnum in residues} 

    atomnames_query = " or ".join([f"name == '{atom}'" for atom in atoms])
    residues_query = " or ".join([f"residue == {resnum}" for resnum in residues])

    if len(atomnames_query) == 0:
        if len(residues_query) == 0:
            filtered_trajectory = trajectory
        else:
            atom_indices_selection = topology.select(residues_query)
            filtered_trajectory = trajectory.atom_slice(atom_indices_selection)
    else:
        if len(residues_query) == 0:
            atom_indices_selection = topology.select(atomnames_query)
            filtered_trajectory = trajectory.atom_slice(atom_indices_selection)
        else:
            atom_indices_selection = topology.select('(' + atomnames_query + ') and (' + residues_query + ')')
            filtered_trajectory = trajectory.atom_slice(atom_indices_selection)
    print("WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed")

    return filtered_trajectory


def filter_traj_to_pdb(trj_file : str, top_file : str, 
                       pdb : str, residues : str | set = {},
                        atoms : set = {}) -> None:
    """
    Filters an input trajectory to only the specified atoms and residues and outputs to pdb

    Filteres an input trajectory that contains all of the atoms in a trajectory to only
    the desired atoms at the desired residues (eg. the atoms necessary to find the 
    center of geometry of a residue) and writes the output to a specified pdb file.
    If residues or atomnames are empty, all residues or atoms are included respectively.

    Parameters
    ----------
    trj_file : str
        path to file of the concatenated trajectory. Should be resampled to the
        1 in 50 frames sampled trajectories for each replicate.
    top_file : str
        path to file of the topology of the molecule
    pdb : str
        path to the output pdb file
    residues : set or str
        1-indexed residue numbers of residues to keep in the trajectory.
        Accepts smart-indexed str representing a list of residues (e.g '1-5,6,39-48')
        If Empty, include all residues.
    atomnames : set 
        atomnames to keep in the trajectory

    Returns
    -------
    None

    See Also
    --------
    filter_traj : Filters an input trajectory to only the specified atoms and residues
    
    Notes
    -----
    Inputed trajectory should have 1-indexed Residue Indices, 
    Outputed trajectory object will be 0-indexed.

    """
    residues = SmartIndexingAction.parse_smart_index(residues)

    filtered_trajectory = filter_traj(trj_file, top_file, residues, atoms)
    filtered_trajectory.save_pdb(pdb)
    print("WARNING: Output file atom, residue, and chain indices are zero-indexed")
    print("Filtered trajectory written to: ", pdb)


def file_convert(trj_file: str, top_file: str, outfile: str) -> None:
    """
    Converts a trajectory input file to a new output type.

    The output file type is determined by the ``outfile`` extension. Uses ``mdtraj.save()`` commands to convert 
    trajectory files to various file types such as ``mdtraj.save_mdcrd()``, ``mdtraj.save_pdb()``, ``mdtraj.save_xyz()``, etc.

    Parameters
    ----------
    trj_file : str
        Path to the file of the concatenated trajectory (eg. .mdcrd file). 
    top_file : str
        Path to the file of the topology of the molecule (.prmtop file).
    outfile : str
        Output filename (include .mdcrd, .pdb, etc.).

    Returns
    -------
    None

    Examples
    --------
    >>> import stacker as st
    >>> import mdtraj as md
    >>> st.file_convert('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                 'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                 'stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.xyz')
    WARNING: Output file atom, residue, and chain indices are zero-indexed
    Trajectory written to: stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.xyz
    >>> md.load_xyz('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.xyz', 
    ...             top='stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop')
    <mdtraj.Trajectory with 10 frames, 12089 atoms, 494 residues, without unitcells at 0x10bb75cd0>

    Notes
    -----
    Output filetype determined from file extension of `output_file` parameter.
    
    See Also
    --------
    mdtraj.load : Load trajectory files
    mdtraj.save : Save md.Trajectory to file
    mdtraj.load_xyz : Load a .xyz trajectory file

    """
    print("WARNING: Output file atom, residue, and chain indices are zero-indexed")
    trajectory = md.load(trj_file, top = top_file)
    trajectory.save(outfile)
    print("Trajectory written to: ", outfile)

class SmartIndexingAction(argparse.Action):
    '''
    Custom argparse action to handle smart indexing of frame numbers.

    Parses a comma-separated list of frame numbers with optional ranges (e.g., '1-20, 34, 25, 50-100')
    and generates a list of individual frame numbers. Modifies the namespace by setting the attribute specified by the 'dest' parameter to the
    list of individual frame numbers.

    Parameters
    ----------
    parser: argparse.ArgumentParser
        The argparse parser object.
    namespace: argparse.Namespace
        The argparse namespace.
    values: str
        The input string containing frame numbers and ranges.
    option_string: str, default=None
        The option string.

    Attributes
    ----------
    dest : str
        The attribute name in the namespace where the parsed list will be stored.

    Methods
    -------
    __call__(parser, namespace, values, option_string=None)
        Parses the provided string of values into a sorted list of integers and
        sets it as an attribute in the namespace.

    Examples
    --------
    >>> parser = argparse.ArgumentParser()
    >>> parser.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze", required=False, action=SmartIndexingAction)
    >>> args = parser.parse_args(["-fl", "1-20,34,25,50-100"])
    >>> print(args.frame_list)
    [1, 2, ..., 20, 34, 25, 50, 51, ..., 100]
    
    '''
    def __call__(self, parser, namespace, values, option_string=None):
        frame_list = []
        for item in values.split(','):
            if '-' in item:
                start, end = map(int, item.split('-'))
                frame_list.extend(range(start, end + 1))
            else:
                frame_list.append(int(item))
        frame_list.sort()
        setattr(namespace, self.dest, frame_list)
    
    @staticmethod
    def parse_smart_index(value : str | set | typing.ArrayLike):
        """
        Checks that an inputted variable is a list that can be smart indexed
        and indexes it if necessary.

        Parameters
        ----------
        value : {str, set, list}
            The input string containing ranges (e.g., "1-5,10,15-20")
            or a set of integers.

        Returns
        -------
        set
            A set of integers parsed from the input.

        Examples
        --------
        >>> import stacker as st
        >>> st.SmartIndexingAction.parse_smart_index('1-5,8,12-17')
        {1, 2, 3, 4, 5, 8, 12, 13, 14, 15, 16, 17}

        Raises
        ------
        ValueError
            If the input is not a string or set.
        """
        if isinstance(value, str):
            parsed_set = set()
            if value == '':
                return {}
            
            for item in value.split(','):
                if '-' in item:
                    start, end = map(int, item.split('-'))
                    parsed_set.update(range(start, end + 1))
                else:
                    parsed_set.add(int(item))
            return parsed_set
        elif isinstance(value, set) or isinstance(value, set) or isinstance(value, dict):
            return value
        elif isinstance(value, list):
            return value
        else:
            raise ValueError("Input must be a string, list, or set of integers.")


if __name__ == "__main__":
    # filter_traj tests
    print('Known Res: 426 = G and 427 = C')
    filtered_traj = filter_traj('testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', {426,427}, {'C2','C4','C6'})
    table, bonds = filtered_traj.topology.to_dataframe()
    print(table)

    ### No Filtering
    print("No Filtering, known trj has 12089 atoms")
    filtered_traj = filter_traj('testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', residues={}, atoms={})
    table, bonds = filtered_traj.topology.to_dataframe()
    print(table)