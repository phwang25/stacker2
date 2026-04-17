"""The Python Functions that control Command Line Options

This module contains the Python routines called when StACKER
is run in the command line.
"""

from .residue_movement import *
from .visualization import *
from .vector import *
from .pairwise_distance import *
from .file_manipulation import *
from .kmeans import *
import argparse
import sys, os
import numpy as np
import pandas as pd
import random
import math

def run_python_command() -> None:
    '''Reads the user's passed in command line and runs the command

    Reads the command line input, runs the associated command with the
    added flags.
    '''
    parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.RawTextHelpFormatter, 
                                     description="Wrapper to run stacker subroutines using the -s flag.\n" + \
                                        "More info on each routine given by `stacker -s ROUTINE -h`")
    global args;
    args, remaining_args = parser.parse_known_args()
    
    # If no flags specified at all
    if not any(vars(args)) and not remaining_args:
        print('usage: stacker -s ROUTINE [-h]\n\n' + \
            'Wrapper to run stacker subroutines using the -s flag.\n' + \
            'More info on each routine given by `stacker -s ROUTINE -h` or `stacker -s ROUTINE --help`\n\n' + \
            'options:\n' +\
            '-s ROUTINE, --script ROUTINE\n' +\
            '            Name of command to use. Options for ROUTINE:\n\n' + \
            '              filter_traj:\n' +\
            '                    filters trajectory and topology files to desired residue numbers and atom names\n' + \
            '              bottaro OR pairwise OR psf:\n' +\
            '                    Create a Pairwise Stacking Fingerprint (PSF): polar plots like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)\n' + \
            '              res_distance:\n' + \
            '                    Get the distance between two residues in a given frame\n' + \
            '              system OR ssf:\n' + \
            '                    Create a System Stacking Fingerprint (SSF) averaged across specified frames\n' + \
            '              stack_events:\n' + \
            '                    Get list of residues with most stacking events (distance closest to 3.5Å)\n' + \
            '              compare:\n' +\
            '                    Get the most changed stacking events between two fingerprints using the outputs of stacker -s stack_events\n\n' +\
            '-h, --help            show this help message and exit\n')
        return

    # help when no script specified
    if ('-s' not in remaining_args and '--script' not in remaining_args and ('--help' in remaining_args or '-h' in remaining_args)): 
        parser.add_argument("-s", "--script", metavar="ROUTINE", help='Name of command to use. Options for ROUTINE:\n\n' + \
                            "  filter_traj:\n\tfilters trajectory and topology files to desired residue numbers and atom names\n" + \
                            "  bottaro OR pairwise OR psf:\n\tCreate Polar Stacking Fingerprint like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)\n" + \
                            "  res_distance:\n\tGet the distance between two residues in a given frame\n" +\
                            "  system OR ssf:\n\tCreate a System Stacking Fingerprint of distances by residue\n" + \
                            "  stack_events:\n\tGet list of residues with most stacking events (distance closest to 3.5Å)\n" +\
                            "  compare:\n\tGet the most changed stacking events between two fingerprints using the outputs of stacker -s stack_events\n",
                             required=True, default='', choices=['filter_traj', 'bottaro', 'pairwise', 'psf', 'res_distance', 'system', 'ssf', 'stack_events', 'compare'])
        parser.add_argument("-h", "--help", help="show this help message and exit", action='help')
        args = parser.parse_args()

    parser.add_argument("-s", "--script", metavar="ROUTINE", help='Name of command to use. Options for ROUTINE:\n\n' + \
                            "  filter_traj:\n\tfilters trajectory and topology files to desired residue numbers and atom names\n" + \
                            "  bottaro OR pairwise OR psf:\n\tCreate Polar Stacking Fingerprint like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)\n" + \
                            "  res_distance:\n\tGet the distance between two residues in a given frame\n" +\
                            "  system OR ssf:\n\tCreate a System Stacking Fingerprint of distances by residue\n" + \
                            "  stack_events:\n\tGet list of residues with most stacking events (distance closest to 3.5Å)\n" +\
                            "  compare:\n\tGet the most changed stacking events between two fingerprints using the outputs of stacker -s stack_events\n",
                             required=True, default='', choices=['filter_traj', 'bottaro', 'pairwise', 'psf', 'res_distance', 'system', 'ssf', 'stack_events', 'compare','test'])
      
    args, remaining_args = parser.parse_known_args()

    # Organizes all possible arguments
    ## Determines which script/subroutine is to be run
    if args.script == 'filter_traj':
        parser.description = 'Filters trajectory and topology files to desired residue numbers and atom names and outputs to a PDB\n\nExamples:\n' +\
                            '[user]$ stacker -s filter_traj -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -o testing/command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -r 426,427 -a C2,C4,C6'
        
        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        required_group.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        required_group.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Filepath of PDB to output to", required=True)
        
        # optional arguments
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=False, action = SmartIndexingAction)
        parser.add_argument("-a", "--atom_names", metavar="ATOM_NAMES", help="Comma-separated list of atom names to filter", required=False, default="C2,C4,C6")

    if args.script == 'bottaro' or args.script == 'pairwise' or args.script == 'psf':
        parser.description = 'Create Polar Stacking Fingerprint (PSF) of the movement of a "viewed residue" from the perspective of a "perspective residue"\nlike those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972). Creates CSV of these values' + \
                                '\n\nExamples:\n' +\
                                '\n[user]$ stacker -s bottaro -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt scatter\n' +\
                            '\n[user]$ stacker -s bottaro -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt heat'
        
        required_group = parser.add_argument_group('Required Arguments')
        parser.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation, if empty then 2-residue PDB expected", required=False, default = '')
        parser.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation, if empty then 2-residue PDB expected", required=False, default = '')
        parser.add_argument("-pdb", "--pdb_input", metavar="PDB_INPUT", help="If trajectory provided: filepath to intermediary PDB file containing two residues, the perspective and viewed nucleotide.\nIf no trajectory given, PDB is expected to already be 2-residue (use -s filter_traj if needed).\nIf empty, will use the same prefix as the trajectory file", required=False, default = '')
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Filepath to output Bottaro values (frame, r, rho, theta) to. If empty, will use the same prefix as the trajectory file.", required=False, default = '')
        required_group.add_argument("-p", "--pers_res", metavar="PERSPECTIVE_RES", help="residue index of the perspective residue whose plane to project onto. 0-/1-indexed based on -i flag (default: 1-indexed)", required=True)
        required_group.add_argument("-v", "--view_res", metavar="VIEWED_RES", help="residue index of the viewed residue whose midpoint will be projected onto perspective residue plane. 0-/1-indexed based on -i flag (default: 1-indexed)", required=True)
        parser.add_argument("-pa", "--pers_atoms", metavar="PERSPECTIVE_ATOMS", help="Comma-separated list of atomnames to use from residue 1 to find center of geometry for perspective nucleotide", required=False, default="C2,C4,C6")
        parser.add_argument("-va", "--view_atoms", metavar="VIEWED_ATOMS", help="Comma-separated list of atomnames to use from residue 2 to find center of geometry for viewed nucleotide", required=False, default="C2,C4,C6")
        parser.add_argument("-i", "--index", metavar="INDEX", type=int, help="index (0-index or 1-index) for perspective/viewed residue numbers (default: 1-indexed)", required=False, default = 1)
        parser.add_argument("-pt", "--plot_type", metavar="PLOT_TYPE", choices = ['scatter', 'heat', ''], help="plot type (scatter or heat) to visualize Bottaro values. If empty string, then just write to csv with no visualization", required=False, default = '')
        parser.add_argument("-po", "--plot_outfile", metavar="PLOT_OUTFILE", help="filename to output plot png to. If empty string, outputs to standard Python vis", required=False, default = '')
        parser.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", default='', help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze,\ngets average distance between residues across these frames\nif empty all frames are used, cannot be used with -f", required=False, action=SmartIndexingAction)
        parser.add_argument("-n", "--no_inter", help="Delete intermediate files after command runs", action = 'store_true', default=False)

    if args.script == 'res_distance':
        parser.description = 'Get the distance between two residues in a given frame\n\n' + \
                                'Examples:\n' +\
                                '[user]$ stacker -s res_distance -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -f 2 --residues 426,427 --atom_names C2,C4,C6'
        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        required_group.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze", required=False)
        required_group.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues, must provide only 2 residues, accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=True, action = SmartIndexingAction)
        parser.add_argument("-b", "--bootstrap", metavar="N_FRAMES", help="Run bootstrap analysis on this residue pairing, sampling N_FRAMES with replacement", required=False, type = int)
        parser.add_argument("-a", "--atom_names", metavar="ATOM_NAMES", help="Comma-separated list of atom names. Three required to get center of geometry for a residue. default = C2,C4,C6", required=False, default="C2,C4,C6")

    if args.script == 'system' or args.script == 'ssf':
        parser.description = 'Creates a System Stacking Fingerprint of the average structure across the chosen frames of a trajectory.' + \
                                '\n\nExamples:\n' +\
                                '\n[user]$ stacker -s system -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -fl 1-2\n' +\
                                '\n[user]$ stacker -s ssf -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -fl 1-2 -g 10 -o testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_nowat_pairwise_avg_1to2.png -d testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_data_1to2.txt\n'
        
        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        required_group.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=False, action = SmartIndexingAction)
        parser.add_argument("-i", "--input", metavar="INPUT_FILE", help="Input .txt file containing per-frame stacking information, in lieu of running stacking fingerprint analysis again.\nTXT file can be created by running `stacker -s system -d OUTPUT_FILE`\n-r flag must match the residues used to create the TXT file")
        frame_group = parser.add_mutually_exclusive_group()
        frame_group.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze, cannot be used with -fl", required=False)
        frame_group.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", default='', help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze,\ngets average distance between residues across these frames\nif empty all frames are used, cannot be used with -f", required=False, action=SmartIndexingAction)
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Filename of output PNG to write plot to. If empty, will output displays to Python visual", default = '', required=False)
        parser.add_argument("-g", "--get_stacking", metavar="N_EVENTS", help="Get list of N_EVENTS residues with most stacking events (distance closest to 3.5Å) in the average structure across all frames.\nPrint to standard output. Equivalent to -s stack_events -n N_EVENTS", type = int, required=False, default = -1)
        parser.add_argument("-d", "--data_output", metavar="OUTPUT_FILE", help="Output the calculated per-frame numpy arrays that create the stacking fingerprint matrix to a file", default = '', required=False)
        parser.add_argument("-B", "--input_B", metavar="INPUT_FILE", help="Input .txt file containing per-frame stacking information for a second fingerprint, creates fingerprint where top left is initial input, bottom right is second fingerprint.\n Used in lieu of running stacking fingerprint analysis again.\nTXT file can be created by running `stacker -s system -d OUTPUT_FILE`\n-r flag must match the residues used to create the TXT file")
        parser.add_argument("-l", "--limits", metavar="LIMITS", help="limits of the color scale, default = (0,7)", default = '0,7')
        parser.add_argument("-y", "--scale_style", metavar="SCALE_STYLE", help=" style of color scale. {bellcurve, gradient}", default = 'bellcurve')
        parser.add_argument("-t", "--threads", metavar="N_THREADS", help="Use multithreading with INT worker threads", type = int, required = False, default = 1)

    if args.script == 'stack_events':
        parser.description = 'Get list of residues with most stacking events (distance closest to 3.5Å) in the stacking fingerprint of the average structure across all frames of a trajectory' + \
                                '\n\nExamples:\n' +\
                                '\n[user]$ stacker -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -f 1 -n 5\n' +\
                                '[user]$ stacker -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-100 -fl 1-10 -n 5\n' +\
                                '[user]$ stacker -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -n 5 -i testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_data_1to2.txt\n' 

        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        required_group.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues to subset trajectory, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=False, action = SmartIndexingAction)
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Output tab-separated txt file to write top stacking events to. If empty, will output displays to standard output", default = '', required=False)
        parser.add_argument("-n", "--n_events", type = int, metavar="N_EVENTS", help="Number of stacking events to display. If -1 display all events", default = -1, required=False)
        parser.add_argument("-i", "--input", metavar="INPUT_FILE", help="Input .txt file containing per-frame stacking information, in lieu of running stacking fingerprint analysis again.\nTXT file can be created by running `stacker -s system -d OUTPUT_FILE`\n-r flag must match the residues used to create the TXT file")
        parser.add_argument("-j", "--include_adjacent", help="Boolean whether to include adjacent residues in the printed output", action = 'store_true', default=False)
        frame_group = parser.add_mutually_exclusive_group()
        frame_group.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze, cannot be used with -fl", required=False)
        frame_group.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", default='', help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze,\ngets average distance between residues across these frames\nif empty all frames are used, cannot be used with -f", required=False, action=SmartIndexingAction)
        parser.add_argument("-t", "--threads", metavar="N_THREADS", help="Use multithreading with INT worker threads", type = int, required = False, default = 1)

    if args.script == 'compare':
        parser.description = 'Print the most changed stacking events between two fingerprints using the outputs of stacker -s stack_events' +\
                                '\n\nExamples:\n' +\
                                '[user]$ stacker -s compare -A /home66/esakkas/STACKER/SCRIPTS/slurmLogs_fingerprint/out_fingerprint_2418986 -B /home66/esakkas/STACKER/SCRIPTS/slurmLogs_fingerprint/out_fingerprint_2418997 -SA _tUAG_aCUA_+1GCU -SB _tUAG_aCUA_+1CGU\n'

        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-A", "--file_A", metavar="FILENAME_A", help = "Filepath to the output log of stacker -s stack_events for the first stacking fingerprint", required = True)
        required_group.add_argument('-B', '--file_B', metavar="FILENAME_B", help = 'Filepath to the output log of stacker -s stack_events for the second stacking fingerprint', required = True)
        required_group.add_argument('-SA', '--source_A', metavar="SOURCE_A", help = 'String describing source of file A, e.g. `_tUAG_aCUA_+1GCU`', required = True)
        required_group.add_argument('-SB', '--source_B', metavar="SOURCE_B", help = 'String describing source of file B, e.g. `_tUAG_aCUA_+1CGU`', required = True)
        parser.add_argument("--nt_cut",type=float,default=4.5,help="COG distance cutoff (Å) for nucleotide–nucleotide pairs (default: 4.5).")
        parser.add_argument("--aa_cut",type=float,default=5.0,help="COG distance cutoff (Å) for pairs involving amino acids (default: 5.0).")

    # help for specific scripts
    if '--help' in remaining_args or '-h' in remaining_args:
        parser.add_argument("-h", "--help", help="show this help message and exit", action='help')
        args = parser.parse_args()

    args = parser.parse_args()
    convert_to_python_command()

def convert_to_python_command() -> None:
    '''Converts a parsed command to use to the correct subroutine and runs the routine

    Converts the specified script to a python command and runs it with the associated inputs
    based on the flags.
    '''
    command = args.script

    if command == 'filter_traj':
        filter_traj_routine()
    elif command == 'bottaro' or command == 'pairwise' or command == 'psf':
        bottaro_routine()
    elif command == 'res_distance':
        res_distance_routine()
    elif command == 'system' or command == 'ssf':
        system_routine()
    elif command == 'stack_events':
        stack_events_routine()
    elif command == 'compare':
        compare_routine()
    else:
        raise InvalidRoutine(args.script + " is not a valid routine")
    
def filter_traj_routine() -> None:
    """
    Executes the routine for filtering an input trajectory file and converting it to a PDB file.

    This function uses provided command-line arguments to call `filter_traj_to_pdb()` 
    with specified inputs, such as the trajectory file, topology file, output file location, 
    residues to retain, and atom names to retain.

    Raises
    ------
    ResEmpty
        If the `--residues` argument is not provided.
    AtomEmpty
        If the `--atom_names` argument is not provided.

    Examples
    --------
    Command-line usage with sample arguments::

        $ stacker 
            -s filter_traj 
            -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd 
            -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop 
            -o testing/command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb 
            -r 426,427 -a C2,C4,C6

    Notes
    -----
    - The `args.residues` argument should be a list of residue indices to retain.
    - The `args.atom_names` argument should be a comma-separated string of atom names to retain.
    - Output directories are created automatically if they do not exist.

    See Also
    --------
    filter_traj_to_pdb : The core function that performs the trajectory filtering and PDB generation.

    """
    if args.residues:
        residues_desired = set(args.residues)
    else:
        raise ResEmpty("Must include a list of residues to keep in the trajectory")

    if args.atom_names:
        atomnames_desired = {atom.strip() for atom in args.atom_names.split(",")}
    else:
        raise AtomEmpty("Must include a list of atom names to keep in the trajectory")

    create_parent_directories(args.output)
    filter_traj_to_pdb(trj_file=args.trajectory, top_file=args.topology, pdb=args.output, residues=residues_desired, atoms=atomnames_desired)

def bottaro_routine() -> None:
    """
    Executes the residue movement routine to generate Bottaro values for a trajectory.

    This routine calculates the `r`, `rho`, and `theta` values for each frame of a PDB trajectory 
    between two specified residues and stores them in a CSV file. Optionally, it visualizes the 
    movement using heatmaps or scatter plots.

    Raises
    ------
    AtomEmpty
        If required atom names or residue indices are not provided.

    Parameters
    ----------
    None
        The routine relies on global `args`, which must be set via command-line arguments or 
        equivalent argument parsing.

    Notes
    -----
    - `args.trajectory` and `args.topology` are required for generating an intermediate PDB if one 
      is not already provided.
    - Supports heatmap and scatter plot visualizations for residue movement.
    - Removes intermediate files (`pdb` and `csv`) if `--no_inter` flag is specified.

    Examples
    --------
    Command-line usage::

        $ stacker -s bottaro 
            -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd 
            -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop 
            -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb 
            -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot.csv 
            -p 426 
            -v 427 
            -pa C2,C4,C6 
            -va C2,C4,C6 
            -pt scatter

    
        $ stacker -s bottaro 
            -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb 
            -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv 
            -p 426 
            -v 427 
            -pa C2,C4,C6 
            -va C2,C4,C6 
            -pt heat

    See Also
    --------
    filter_traj_to_pdb : Filters trajectory and generates an intermediate PDB.
    write_bottaro_to_csv : Writes residue movement calculations to a CSV.
    visualize_two_residue_movement_heatmap : Creates a heatmap for residue movement.
    visualize_two_residue_movement_scatterplot : Creates a scatter plot for residue movement.


    """
    trj_prefix = args.trajectory.rsplit('.', 1)[0] 

    if args.pdb_input == '':
        pdb_filename = trj_prefix + '.pdb'
    else:
        pdb_filename = args.pdb_input

    if args.output == '':
        prefix = pdb_filename.rsplit('.', 1)[0] 
        output_name = prefix + '.csv'
    else:
        output_name = args.output

    if args.pers_atoms:
        perspective_atom_names = {res.strip() for res in args.pers_atoms.split(",")}
    else:
        raise AtomEmpty("Must include a list of atom names to define Pespective Residue center of geometry")
    
    if args.view_atoms:
        viewed_atom_names = {res.strip() for res in args.view_atoms.split(",")}
    else:
        raise AtomEmpty("Must include a list of atom names to define Viewed Residue center of geometry")

    if args.pers_res is not None:
        pers_res_num = int(args.pers_res)
    else:
        raise AtomEmpty("Must include a 1-indexed residue index for the perspective residue")
    
    if args.view_res is not None:
        view_res_num = int(args.view_res)
    else:
        raise AtomEmpty("Must include a 1-indexed residue index for the perspective residue")

    if args.frame_list:
        frame_list = set(args.frame_list)
    else:
        frame_list = {}

    create_parent_directories(pdb_filename)
    if args.trajectory and args.topology:
        filter_traj_to_pdb(trj_file=args.trajectory, top_file=args.topology, pdb=pdb_filename,
                           residues={pers_res_num,view_res_num}, atoms=perspective_atom_names.union(viewed_atom_names))
    
    create_parent_directories(output_name)
    
    write_bottaro_to_csv2(pdb=pdb_filename, 
                         outcsv=output_name, pers_res=pers_res_num, view_res=view_res_num,
                         res1_atoms=tuple(perspective_atom_names), 
                         res2_atoms=tuple(viewed_atom_names), index = args.index)

    # write_bottaro_to_csv(pdb=pdb_filename, 
    #                      outcsv=output_name, pers_res=pers_res_num, view_res=view_res_num,
    #                      res1_atoms=tuple(perspective_atom_names), 
    #                      res2_atoms=tuple(viewed_atom_names), index = args.index)
    
    # if args.plot_type == 'heat':
    #     create_parent_directories(args.plot_outfile)
    #     visualize_two_residue_movement_heatmap(output_name, plot_outfile=args.plot_outfile, frame_list = frame_list)
    # elif args.plot_type == 'scatter':
    #     create_parent_directories(args.plot_outfile)
    #     visualize_two_residue_movement_scatterplot(output_name, plot_outfile=args.plot_outfile, frame_list = frame_list)

    # if args.no_inter:
    #     if args.trajectory and args.topology: #intermediate pdb created
    #         os.remove(pdb_filename)
    #     if os.path.exists(output_name):
    #         os.remove(output_name)

def res_distance_routine() -> None:
    """
    Calculates the distance between the centers of geometries of two specified residues.

    This routine filters the trajectory based on the specified residues and atom names, 
    and calculates the distances for either a single frame or multiple frames if bootstrapping is enabled.

    Raises
    ------
    ResEmpty
        If fewer or more than two residues are specified, or if no residues are provided.
    AtomEmpty
        If no atom names are provided.

    Parameters
    ----------
    None
        The routine relies on global `args` for all input values.

    Notes
    -----
    - The `args.residues` argument must contain exactly two residues.
    - The `args.atom_names` argument is required to specify the atoms involved in the calculation.
    - If the `--bootstrap` argument is provided, distances are calculated across a specified number of randomly sampled frames.
    - If no bootstrap is performed, the distance is calculated for a single frame.

    Examples
    --------
    Command-line usage::

        $ stacker -s res_distance 
            -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd 
            -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop 
            -f 2 
            --residues 426,427 
            --atom_names C2,C4,C6

    Outputs
    -------
    If bootstrapping is enabled:
        - Prints the mean distance, standard deviation, and percentiles for the distances.
    Otherwise:
        - Prints the calculated distance for the specified frame.

    See Also
    --------
    filter_traj : Filters the trajectory for the specified residues and atom names.
    calculate_residue_distance : Calculates the distance between residues for a given frame.
    
    """
    if len(args.residues) != 2:
        raise ResEmpty("Must include only 2 residues")
    elif args.residues:
        residues_desired = set(args.residues)
    else:
        raise ResEmpty("Must include a list of residues to keep in the trajectory")

    if args.atom_names is not None:
        atomnames_desired = {atom.strip() for atom in args.atom_names.split(",")}
    else:
        raise AtomEmpty("Must include a list of atom names to keep in the trajectory")

    block_printing()
    filtered_trj = filter_traj(trj_file=args.trajectory, top_file=args.topology, residues=residues_desired, atoms=atomnames_desired)
    if args.bootstrap:
        n_frames = len(filtered_trj)
        frames = [random.randint(0, n_frames-1) for _ in range(args.bootstrap)]
    else:
        frames = [args.frame-1]

    i = 0
    residues_desired = list(residues_desired)
    res_distances = []
    for frame in frames:
        trj_frame = filtered_trj[frame]
        # Correct that calculate_residue_distance res_nums are 1-indexed
        print(i)
        i+=1
        distance_vector = calculate_residue_distance(trj=trj_frame, res1=int(residues_desired[0]), res2=int(residues_desired[1]), res1_atoms=tuple(atomnames_desired), res2_atoms=tuple(atomnames_desired))
        enable_printing()
        res_distances.append(distance_vector.magnitude())

    if args.bootstrap:
        print("Bootstrap Mean Distance:", np.mean(res_distances))
        print("Bootstrap Standard Deviation:", np.std(res_distances))
        print("Bootstrap Percentile:", np.percentile(res_distances, [2.5,5.0, 95.0, 97.5]))
    else:
        print(res_distances[0])

def system_routine() -> None:
    """
    Runs the System Stacking Fingerprint (SSF) routine to generate a single SSF for a specified frame or range of frames.

    This routine processes a trajectory to create SSFs based on inter-residue distances and outputs either visualizations 
    or raw data for further analysis. The user can provide input fingerprints or generate them from trajectory files.

    Raises
    ------
    FileNotFoundError
        If the specified trajectory or topology file does not exist.
    ValueError
        If the provided input files have incompatible dimensions or unexpected formatting.

    Parameters
    ----------
    None
        The routine relies on global `args` for all input values.

    Notes
    -----
    - If `args.input` is provided, the routine will use pre-calculated fingerprint data.
    - If `args.residues` is specified, only those residues are considered in the calculation.
    - SSFs can be created for individual frames, a list of frames, or all frames in the trajectory.
    - If `args.input_B` is provided, two SSFs are combined for analysis.

    Outputs
    -------
    - Saves raw SSF data to `args.data_output` if specified.
    - Generates visualizations of the SSF data as videos saved to `args.output`.
    - Outputs information about the most stacked residues if `args.get_stacking` is specified.

    Examples
    --------
    Command-line usage::

        $ stacker -s system 
            -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd 
            -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop 
            -r 90-215 
            -fl 1-2 

        $ stacker -s ssf 
            -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd 
            -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop 
            -r 90-215 
            -fl 1-2 
            -g 10 
            -o testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_nowat_pairwise_avg_1to2.png 
            -d testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_data_1to2.txt.gz

    See Also
    --------
    filter_traj : Filters the trajectory for the specified residues.
    get_residue_distance_for_trajectory : Computes inter-residue distances for a trajectory.
    combine_frames : Combines two SSFs into a single array for comparison.
    display_arrays_as_video : Visualizes SSF arrays and saves them as video files.
    
    """
    if args.residues:
        residues_desired = set(args.residues)
    else:
        residues_desired = {}

    if args.frame_list:
        frame_list = args.frame_list
    else:
        frame_list = []

    # Reduce the trajectory for specific residues (saves lot of time!)
    trj_sub = filter_traj(trj_file=args.trajectory, top_file=args.topology, residues=residues_desired)

    if args.input: # For loading a previously run fingerprint
        print("Loaded fingerprint data from:", args.input)
        frames = np.load(args.input)
    
    elif args.frame_list: # For ssf of desired frames
        frames = get_residue_distance_for_trajectory(trj_sub, frame_list, threads = args.threads)
    
    elif args.frame: # For single frame ssf
        frames = np.array([get_residue_distance_for_frame(trj_sub, args.frame)])
    
    else: # For trajectory ssf 
        frames = get_residue_distance_for_trajectory(trj_sub, [i for i in range(1,trj_sub.n_frames+1)], threads = args.threads)

    if args.data_output:
        print(f"Number of SSFs made:{frames.shape}")
        np.save(args.data_output, frames)

    if args.limits:
        scale_limits = tuple(float(i) for i in args.limits.replace('(','').replace(')','').replace('...', '').split(','))

    # Average the frames
    avg_frames = [get_frame_average(frames)]

    if args.input_B:
        print("Loaded second fingerprint data from:", args.input_B)
        frames_B = np.load(args.input_B)
        avg_frames_B = [get_frame_average(frames_B)]
        avg_frames = [combine_frames(avg_frames[0], avg_frames_B[0])]
        print(avg_frames)

    if args.get_stacking:
        get_top_stacking(trj_sub, avg_frames[0], csv = '', n_events = args.get_stacking)

    sorted_res = list(residues_desired)
    sorted_res.sort()

    create_parent_directories(args.output)
    
    # Save the fingerprint image
    display_arrays_as_video(avg_frames, sorted_res, seconds_per_frame=10, outfile=args.output, scale_limits=scale_limits, scale_style=args.scale_style)

def combine_frames(frames_A, frames_B):
    """
    Combines two 2D numpy arrays (frames_A and frames_B) into a single array.

    This function takes two 2D numpy arrays of the same shape and combines them
    into a new array. The upper triangular part (excluding the diagonal) of the 
    resulting array is filled with the corresponding elements from frames_A, 
    while the lower triangular part (including the diagonal) is filled with the 
    corresponding elements from frames_B.

    Parameters
    ----------
    frames_A : numpy.ndarray
        A 2D numpy array.
    frames_B : numpy.ndarray 
        A 2D numpy array of the same shape as frames_A.

    Returns
    -------
    numpy.ndarray 
        A new 2D numpy array with combined elements from frames_A and frames_B.

    Examples
    --------
    >>> import stacker as st
    >>> import numpy as np
    >>> frames_A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> frames_B = np.array([[9, 8, 7], [6, 5, 4], [3, 2, 1]])
    >>> st.combine_frames(frames_A, frames_B)
    array([[9., 2., 3.],
           [6., 5., 6.],

           [3., 2., 1.]])

    """
    Am, An = frames_A.shape
    array_to_fill = np.zeros((Am,An))
    for i in range(Am):
        array_to_fill[i,i:] = frames_A[i,i:] 
    for j in range(An):
        array_to_fill[j:,j] = frames_B[j:,j]
    return array_to_fill

def stack_events_routine() -> None:
    """
    Identifies residue pairs with the most π-stacking interactions.

    This routine analyzes a molecular dynamics trajectory and identifies residue pairs 
    with the smallest center of geometry (COG) distances, typically indicating strong 
    π-stacking. It outputs the top `n_events` residue pairings based on user input.

    Raises
    ------
    FileNotFoundError
        If the trajectory or topology files are not found.
    ValueError
        If input parameters are invalid or incompatible with the data.

    Parameters
    ----------
    None
        This function relies on global `args` for input values.

    Outputs
    -------
    - CSV file listing the top stacking events (`args.output`), if specified.

    Examples
    --------
    Command-line usage::

        $ stacker -s stack_events 
            -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd 
            -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop 
            -r 90-215 
            -f 1 
            -n 5

        $ stacker -s stack_events 
            -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd 
            -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop 
            -r 90-100 
            -fl 1-10 
            -n 5

        $ stacker -s stack_events 
            -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd 
            -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop 
            -r 90-215 
            -n 5 
            -i testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_data_1to2.txt

    See Also
    --------
    filter_traj : Filters trajectory data based on residue selection.
    get_top_stacking : Extracts top stacking events from a trajectory.

    """
    if args.residues:
        residues_desired = set(args.residues)
    else:
        residues_desired = {}

    trj_sub = filter_traj(trj_file=args.trajectory, top_file=args.topology, residues=residues_desired)

    if args.input:
        loaded_arr = np.loadtxt(args.input)
        frames = loaded_arr.reshape(loaded_arr.shape[0], math.isqrt(loaded_arr.shape[1]), math.isqrt(loaded_arr.shape[1]))
        frame = get_frame_average(frames)
    elif args.frame:
        frame = get_residue_distance_for_frame(trj_sub, frame = args.frame)
    else:
        frames = get_residue_distance_for_trajectory(trj_sub, args.frame_list, threads = args.threads)
        frame = get_frame_average(frames)

    get_top_stacking(trj_sub, frame, csv = args.output, n_events = args.n_events, include_adjacent = args.include_adjacent)

def compare_routine() -> None:
    """
    Compares pi-stacking events between two trajectories.

    Uses two annotated stack_events TSVs (from get_top_stacking) and prints
    residue pairs with the largest change in average COG distance between
    the two systems, applying separate cutoffs for nt-nt and AA-involving pairs.
    """
    file1 = args.file_A
    file2 = args.file_B
    file1_source = args.source_A 
    file2_source = args.source_B

    nt_cut = args.nt_cut
    aa_cut = args.aa_cut
    # nt_cut = 4.5
    # aa_cut = 5.0

    def preprocess_df(df: pd.DataFrame) -> pd.DataFrame:
        """
        Canonicalize residue pairs so that Res1 <= Res2 and keep
        Res1_* / Res2_* columns in sync when swapping.
        """
        df = df.copy()
        swap_mask = df["Res1"] > df["Res2"]

        # swap numeric residue IDs
        df.loc[swap_mask, ["Res1", "Res2"]] = (
            df.loc[swap_mask, ["Res2", "Res1"]].to_numpy()
        )

        # swap per-residue metadata columns if present
        paired_suffixes = ["_name", "_class"]
        for suf in paired_suffixes:
            c1 = "Res1" + suf
            c2 = "Res2" + suf
            if c1 in df.columns and c2 in df.columns:
                tmp = df.loc[swap_mask, c1].copy()
                df.loc[swap_mask, c1] = df.loc[swap_mask, c2]
                df.loc[swap_mask, c2] = tmp

        # Ensure Pair_type exists
        if "Pair_type" not in df.columns:
            # fall back to simple infer from classes if available
            if "Res1_class" in df.columns and "Res2_class" in df.columns:
                def _infer_pair_type(row):
                    c1 = row["Res1_class"]
                    c2 = row["Res2_class"]
                    if c1 == "nt" and c2 == "nt":
                        return "nt-nt"
                    if c1.startswith("aa") and c2.startswith("aa"):
                        return "aa-aa"
                    return "aa-nt"
                df["Pair_type"] = df.apply(_infer_pair_type, axis=1)
            else:
                # last-resort: treat everything as nt-nt
                df["Pair_type"] = "nt-nt"

        return df

    # Load files (annotated stack_events output)
    try:
        data1 = pd.read_csv(file1, sep="\t")
        data2 = pd.read_csv(file2, sep="\t")
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Error reading stack_events files: {e}")

    data1 = preprocess_df(data1)
    data2 = preprocess_df(data2)

    # Merge on residue pairs; suffix other columns by source tags
    merged_data = pd.merge(
        data1,
        data2,
        on=["Res1", "Res2"],
        suffixes=[file1_source, file2_source],
        how="inner",
    )

    # Unify Pair_type from one side (should be identical)
    pt_col1 = "Pair_type" + file1_source
    pt_col2 = "Pair_type" + file2_source

    if pt_col1 in merged_data.columns:
        merged_data["Pair_type"] = merged_data[pt_col1]
    elif "Pair_type" in merged_data.columns:
        # unsuffixed case
        pass
    else:
        merged_data["Pair_type"] = "nt-nt"

    # Compute absolute discrepancy in Avg_Dist
    col1 = "Avg_Dist" + file1_source
    col2 = "Avg_Dist" + file2_source
    if col1 not in merged_data.columns or col2 not in merged_data.columns:
        raise ValueError(
            f"Expected columns '{col1}' and '{col2}' in merged data; "
            f"found: {list(merged_data.columns)}"
        )

    merged_data["Discrepancy"] = (merged_data[col1] - merged_data[col2]).abs()

    # Apply separate distance cutoffs for nt-nt vs AA-involving pairs
    is_nt_nt = merged_data["Pair_type"] == "nt-nt"
    is_aa_involved = ~is_nt_nt

    subset_data = merged_data[
        (is_nt_nt & ((merged_data[col1] < nt_cut) | (merged_data[col2] < nt_cut)))
        |
        (is_aa_involved & ((merged_data[col1] < aa_cut) | (merged_data[col2] < aa_cut)))
    ].copy()

    subset_data = subset_data.sort_values(by="Discrepancy", ascending=False)

    # Print everything to stdout like the old version
    # (You can add an --output option later if you want to write TSV.)
    if subset_data.empty:
        print("No residue pairs passed the distance cutoffs.")
    else:
        print(subset_data.to_string(index=False))

class InvalidRoutine(Exception):
    '''Specified command line routine was invalid

    This is raised when the `stacker -s ROUTINE` was given
    an invalid ROUTINE.

    Examples
    --------
    Command Line::

        $ stacker -s blah

    '''
    pass

class ResEmpty(Exception):
    """No Residues to subset to were provided

    This is raised when a routine is given an invalid
    number of residues to subset to.

    Examples
    --------

    Command Line::

        $ stacker -s filter_traj -r -a C2,C4,C6
    
    See Also
    --------
    fiter_traj_routine
    bottaro_routine
    """
    pass

class AtomEmpty(Exception):
    """No Atomnames to subset to were provided

    This is raised when a routine is given an invalid
    number of atomnames to subset to.

    Examples
    --------

    Command Line::

        $ stacker -s filter_traj -a -r 426,427  
    
    See Also
    --------
    fiter_traj_routine
    bottaro_routine
    """
    pass

class FrameEmpty(Exception):
    """No Frames present in trajectory"""
    pass

def block_printing():
    '''Disable printing to standard output
    
    References
    ----------
    [1] https://stackoverflow.com/a/8391735

    '''
    sys.stdout = open(os.devnull, 'w')

def enable_printing():
    '''Enable printing to standard output
    
    References
    ----------
    [1] https://stackoverflow.com/a/8391735
    '''
    sys.stdout = sys.__stdout__

if __name__ == '__main__':
    run_python_command()