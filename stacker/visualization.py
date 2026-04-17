"""
Visualize the SSFs, PSFs, and other analyses

This module includes the functions used to visualize plots,
including SSFs and PSFs. The data inputs to these plot functions
are provided by the other modules.
"""

import os
import functools
import numpy as np
from numpy import typing
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from seaborn import kdeplot
from .file_manipulation import SmartIndexingAction

# Optional dependency for PSF v2 plotting.
# StACKER already depends on mdtraj for trajectory/PDB handling, but we
# guard the import so documentation builds (or minimal installs) don't crash.
try:
    import mdtraj as md
except Exception:  # pragma: no cover
    md = None


def create_parent_directories(outfile_prefix : str) -> None:
    '''
    Creates necessary parent directories to write an outfile given a prefix
    
    Parameters
    ----------
    outfile_prefix : str
        The filepath of the output file, including the path where the file will be saved.

    Examples
    --------
    >>> create_parent_directories('/path/to/output/file.txt')
    # This will create the directories '/path/to/output' if they do not exist.
    >>> create_parent_directories('output/file.txt')
    # This will create the directory 'output' if it does not exist.
    '''
    dir_name = os.path.dirname(outfile_prefix)
    if dir_name == '': dir_name = '.'
    os.makedirs(dir_name, exist_ok=True)

def create_axis_labels(res_indices, tick_distance: int = 10, min_tick_spacing: int = 3, max_ticks: int = 60):
    """
    Improved axis labeling for SSFs with many disjoint residue blocks.

    - Always labels start/end of each consecutive block
    - Optionally labels internal positions every `tick_distance`
    - Deduplicates + enforces minimum tick spacing
    - Caps number of ticks to avoid label explosions
    """
    res_indices = list(res_indices)
    n = len(res_indices)
    if n < 1:
        raise NoResidues("pairwise analysis must include at least one residue")

    # Identify blocks in *matrix index space*
    blocks = []
    start = 0
    for i in range(1, n):
        if res_indices[i] != res_indices[i - 1] + 1:
            blocks.append((start, i - 1))
            start = i
    blocks.append((start, n - 1))

    # Candidate ticks (matrix indices)
    ticks = set()
    priority = set()  # block boundaries get priority

    for s, e in blocks:
        ticks.add(s); ticks.add(e)
        priority.add(s); priority.add(e)

        # Internal ticks for long blocks
        if tick_distance and (e - s) >= tick_distance:
            j = s + tick_distance
            while j < e:
                ticks.add(j)
                j += tick_distance

    ticks = sorted(ticks)

    # Enforce minimum spacing, but never drop boundary ticks
    kept = []
    for t in ticks:
        if not kept:
            kept.append(t); continue
        if (t - kept[-1]) >= min_tick_spacing or t in priority:
            kept.append(t)

    # Cap number of ticks if still too many (thin uniformly, keep ends)
    if max_ticks and len(kept) > max_ticks:
        step = int(np.ceil(len(kept) / max_ticks))
        kept = [kept[0]] + kept[1:-1:step] + [kept[-1]]

    labels = [res_indices[t] for t in kept]
    return kept, labels

# TODO Old version (Remove later)
'''def create_axis_labels(res_indices: typing.ArrayLike, tick_distance: int = 10) -> list:
    """
    Designates the axis labels to use in the SSF plot.

    Helper function for visualizing SSFs.
    Returns the x-axis tick positions and labels to use on the ticks based on the 
    residues used in a specific SSF analysis. Meant to be used when many 
    disjoint sets of residue indices are used. Ticks will be present every `tick_distance` 
    residues in a collection of adjacent residues, and a tick will exist at both
    ends of any consecutive residue sequence.

    Parameters
    ----------
    res_indices : list
        The list of residue indices used in the pairwise analysis.
        Parameter `residue_desired` passed to `filter_traj()`
    tick_distance : int, default = 10
        Distance between ticks in blocks of consecutive residues.

    Returns
    -------
    tick_locations : list
        List of tick positions (0-based) to place labels on in the axes.
    tick_labels : list
        List of labels to place at the adjacent tick locations.

    See Also
    --------
    filter_traj : Filters an input trajectory to desired residues

    Examples
    --------
    Residues 0-12,98-100 were used. The SSF will label 0,10,12,98,100,
    provided in the second returned list. The first returned list gives 
    the positions on the axes to place each label.

    >>> create_axis_labels([0,1,2,3,4,5,6,7,8,9,10,11,12,98,99,100])
    [0, 10, 12, 13, 15], [0, 10, 12, 98, 100]

    >>> create_axis_labels([94,95,96,97,98,99,100,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428])
    [0, 6, 7, 17, 27], [94, 100, 408, 418, 428]
    """
    n_residues = len(res_indices)

    if n_residues < 1: raise NoResidues("pairwise analysis must include at least one residue")

    tick_locations = [0]
    tick_labels = [res_indices[0]]

    res_sequence_length = 1

    for i in range(1, n_residues):
        if res_indices[i] == res_indices[i-1] + 1:
            res_sequence_length += 1

        if res_indices[i] > res_indices[i-1] + 1:
            tick_locations += [i-1, i]
            tick_labels += [res_indices[i-1], res_indices[i]]
            res_sequence_length = 1
        elif res_sequence_length == tick_distance+1:
            tick_locations += [i]
            tick_labels += [res_indices[i]]
            res_sequence_length = 1
    
    if n_residues-1 not in tick_locations:
        tick_locations += [n_residues-1]
        tick_labels += [res_indices[n_residues-1]]
    
    return tick_locations, tick_labels
'''

def display_arrays_as_video(numpy_arrays : list | typing.ArrayLike, res_indices: typing.ArrayLike | str, 
                            seconds_per_frame: int = 10, tick_distance: int = 10, min_tick_spacing: int = 3, max_ticks: int = 60,
                            outfile_prefix: str = '', scale_limits: tuple = (0, 7), outfile: str = '',
                            scale_style: str = 'bellcurve', xy_line: bool = True, **kwargs) -> None:
    """
    display_arrays_as_video(
        ssfs,
        res_indices,
        seconds_per_frame = 10,
        tick_distance = 10,
        outfile_prefix = '',
        scale_limits = (0,7),
        outfile = '',
        scale_style = 'bellcurve',
        xy_line = True,
        **kwargs
    )

    Displays SSF data to output or writes SSF as a PNG

    Visualizes the data for an SSF for a trajectory or a single frame.
    Takes an SSF array outputted from `get_residue_distance_for_frame`,
    `get_residue_distance_for_trajectory`,
    or `system_stacking_fingerprints` and treats them as frames 
    in a video, filling in a grid at position i, j by the value 
    at i, j in the array.

    Parameters
    ----------
    ssfs : array_like
        List or array of 2D NumPy arrays representing SSFs, 
        output of ``system_stacking_fingerprints``
    res_indices : list or str
        The list of residue indices used in the pairwise analysis.
        Parameter `residue_desired` passed to `filter_traj()`
        Accepts smart-indexed str representing a list of residues (e.g '1-5,6,39-48')
    seconds_per_frame : int, default = 10
        Number of seconds to display each matrix for.
    tick_distance : int, default = 10
        Distance between ticks in blocks of consecutive residues.
    outfile : str
        Image output filepath to write a single SSF to. Format inferred from file extension.
        png, pdf, ps, eps, and svg supported.
    outfile_prefix : str
        Prefix for image filepath to write multiple frames to. Format inferred from file extension.
        png, pdf, ps, eps, and svg supported.
    scale_limits : tuple, default = (0, 7)
        Limits of the color scale.
    scale_style : {'bellcurve', 'gradient'}, default = 'bellcurve'
        Style of color scale. 
    xy_line : bool, default = True
        Draw x = y line to separate matrix halves.
    **kwargs : dict, optional
        Additional keyword arguments to customize the plot:

        - fontsize : int, default = 10
            Font size for all text elements.
        - fig_width : float, default = 8
            Width of the figure in inches.
        - fig_height : float, default = 8
            Height of the figure in inches.
        - title_fontsize : int, default = fontsize
            Font size for the title.
        - legend_fontsize : int, default = fontsize
            Font size for the legend.
        - cb_fontsize : int, default = fontsize
            Font size for the colorbar tick labels.
        - xaxis_fontsize : int, default = fontsize
            Font size for the x-axis tick labels.
        - yaxis_fontsize : int, default = fontsize
            Font size for the y-axis tick labels.

    Returns
    -------
    None
        Displays video of NumPy arrays.

    See Also
    --------
    create_axis_labels : Designates the axis labels to use in the SSF plot.
    get_residue_distance_for_frame : Calculates System Stacking Fingerprint (SSF) between all residues in a given frame.
    get_residue_distance_for_trajectory : get SSF data for all frames of a trajectory
    system_stacking_fingerprints : Alias for this `get_residue_distance_for_trajectory`
    display_ssfs : Alias for this function.
    
    Examples
    --------
    >>> import stacker as st
    >>> import mdtraj as md
    >>> trajectory_file = 'stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    >>> topology_file = 'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
    >>> trj = md.load(trajectory_file, top = topology_file)
    >>> residue_selection_query = 'resi 90 to 215'
    >>> frames_to_include = [1,2,3,4,5]
    >>> trj_sub = trj.atom_slice(trj.top.select(residue_selection_query))
    >>> resSeqs = [res.resSeq for res in trj_sub.topology.residues]
    >>> frames = st.get_residue_distance_for_trajectory(trj_sub, frames_to_include, threads = 5)
    Frame 2 done.
    Frame 3 done.
    Frame 1 done.
    Frame 5 done.
    Frame 4 done.
    >>> st.display_arrays_as_video([st.get_frame_average(frames)], resSeqs, seconds_per_frame=10)
    # Displays SSF for each frame of this trajectory to standard output

    """
    fontsize = kwargs.get('fontsize', 10)
    plt.rcParams.update({'font.size': fontsize})

    orange_colormap = plt.cm.get_cmap('Oranges_r', 100)

    if scale_style == 'gradient':
        newcolors = np.vstack((orange_colormap(np.linspace(1, 1, 1)), orange_colormap(np.linspace(0, 0, 128)),
                        orange_colormap(np.linspace(0, 1, 128))))
        newcmp = plt.cm.colors.ListedColormap(newcolors, name='OrangeBellcurve')
    elif scale_style == 'bellcurve':
        newcolors = np.vstack((orange_colormap(np.linspace(1, 0, 128)), orange_colormap(np.linspace(0, 1, 128))))
        newcmp = plt.cm.colors.ListedColormap(newcolors, name='OrangeBellcurve')
    
    fig = plt.figure(figsize=(kwargs.get('fig_width', 8), kwargs.get('fig_height', 8)))
    ax = fig.add_subplot(111)
    plt.ion()
    frame_num = 1
    for hist in numpy_arrays:
        ax.clear()
        vmin, vmax = scale_limits
        neg = ax.imshow(hist, cmap=newcmp, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_title('Distance Between Residues Center of Geometries', fontsize=kwargs.get('title_fontsize', fontsize))
        ax.set_xlabel('Residue Index')  
        ax.xaxis.set_label_position('top')  
        ax.set_ylabel('Residue Index')  
        colorbar = fig.colorbar(neg, ax=ax, location='right', anchor=(0, 0.3), shrink=0.7)
        colorbar.ax.set_title('Center of\nGeometry\nDist. (Å)', fontsize=kwargs.get('legend_fontsize', fontsize))
        colorbar.ax.tick_params(labelsize=kwargs.get('cb_fontsize', fontsize))


        res_indices = SmartIndexingAction.parse_smart_index(res_indices)
        res_indices = list(res_indices)
        ticks, labels = create_axis_labels(res_indices, tick_distance, min_tick_spacing, max_ticks)
        plt.xticks(ticks, labels, rotation='vertical', fontsize=kwargs.get('xaxis_fontsize', fontsize))
        plt.yticks(ticks, labels, fontsize=kwargs.get('yaxis_fontsize', fontsize))
        ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)


        gap_threshold = kwargs.get("gap_threshold", 10)  # you can pass gap_threshold=10 in your call
        for i in range(1, len(res_indices)):
            gap = res_indices[i] - res_indices[i-1]
            if gap >= gap_threshold:
                b = i - 0.5  # boundary between matrix indices i-1 and i
                ax.axvline(b, linewidth=0.25, alpha=0.12)
                ax.axhline(b, linewidth=0.25, alpha=0.12)


        # separate ticks by region
        # last_res = 1
        # last_label = 0
        # long_tick_region = False
        # for res_i, label_i, tick_object, y_object in zip(ticks, labels, ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()): 
        #     if res_i == last_res+1 and not label_i == last_label+1:
        #         long_tick_region = not long_tick_region

        #     if long_tick_region:
        #         tick_object.set_pad(22.)
        #         tick_object.tick2line.set_markersize(22.)
        #         y_object.set_pad(22.)
        #         y_object.tick1line.set_markersize(22.)
        #     else:
        #         tick_object.set_pad(2.)
        #         tick_object.tick2line.set_markersize(3.)
        #         y_object.set_pad(2.)
        #         y_object.tick1line.set_markersize(3.)
            
        #     last_res = res_i
        #     last_label = label_i

        if xy_line:
            lims = [np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
                    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
                    ]
            ax.plot(lims, lims, 'k-', zorder=0, alpha=0.75, linewidth=0.5) 

        plt.pause(seconds_per_frame)
        if outfile_prefix:
            plt.savefig(outfile_prefix + "frame" + str(frame_num) + ".png")
        elif outfile:
            plt.savefig(outfile)
        colorbar.remove()
        frame_num += 1

@functools.wraps(display_arrays_as_video)
def display_ssfs(*args, **kwargs):
    return display_arrays_as_video(*args, **kwargs)

display_ssfs.__doc__ = f"""
Alias for `display_arrays_as_video()`.

{display_arrays_as_video.__doc__}
"""


# ————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

def set_polar_grid(**kwargs) -> mpl.projections.polar.PolarAxes:
    """
    Set up axes for PSF.

    Creates polar plot background for two-residue movement comparison
    with theta 0 to 360, a radial maximum of 15 Angstroms, and a visualization 
    of the perspective residue at the center.

    Parameters
    ----------
    **kwargs : dict, optional
        Additional keyword arguments to customize the plot:

        - fontsize : int, default = 10
            Font size for all text elements.
        - fig_width : float, default = 5
            Width of the figure in inches.
        - fig_height : float, default = 5
            Height of the figure in inches.
        - title_fontsize : int, default = fontsize
            Font size for the title.
        - legend_fontsize : int, default = fontsize
            Font size for the legend.
        - cb_fontsize : int, default = fontsize
            Font size for the colorbar labels.
        - xaxis_fontsize : int, default = fontsize
            Font size for the x-axis labels.
        - yaxis_fontsize : int, default = fontsize
            Font size for the y-axis labels.

    Returns
    -------
    ax : matplotlib.projections.polar.PolarAxes
        Axis object for the created polar plot.

    """
    fontsize = kwargs.get('fontsize', 10)
    plt.rcParams.update({'font.size': fontsize})

    fig = plt.figure(figsize=(kwargs.get('fig_width', 5), kwargs.get('fig_height', 5)))
    ax = fig.add_subplot(111, polar=True)

    ax.set_xticks(np.pi/180. * np.linspace(0,  360, 3, endpoint=False))
    ax.set_xticklabels([r"$\theta=0^\circ$", r"$\theta=120^\circ$", r"$\theta=240^\circ$"], 
                       fontsize=kwargs.get('xaxis_fontsize', fontsize))
    ax.tick_params(pad=-10)

    ax.set_rlim(0, 10)
    ax.set_rticks(np.linspace(0, 10, 3, endpoint=True))
    ax.tick_params(axis='y', labelsize=kwargs.get('yaxis_fontsize', fontsize))
    ax.set_rlabel_position(180) 
    plt.text(x=np.radians(178), y=12, s=r"$\rho\text{ }(\AA)$", ha="center", va='center', fontsize=kwargs.get('fontsize', fontsize))
    plt.text(x=0, y=2, s="C2", ha="center", va='center', fontsize=kwargs.get('fontsize', fontsize))
    plt.text(x=np.radians(240), y=2, s="C4", ha="center", va='center', fontsize=kwargs.get('fontsize', fontsize))
    plt.text(x=np.radians(120), y=1.8, s="C6", ha="center", va='center', fontsize=kwargs.get('fontsize', fontsize))

    ax.grid(color='gray', linestyle='--', linewidth=0.5)
    return ax

def visualize_two_residue_movement_scatterplot(csv: str, plot_outfile: str = '', frame_list: set = {}, **kwargs) -> None:
    """
    Creates scatterplot of two-residue movement relative to each other.

    Takes the data created in residue_movement and visualizes it as a polar coordinate
    scatterplot similar to the Figure D link in Proposal Feature 4.

    Parameters
    ----------
    csv : str
        Filepath to csv file containing data on the movement
        of two residues relative to each other (r, rho, and theta values). Created
        in residue_movement.
    plot_outfile : str
        Filepath of the image file to write to. Format inferred from file extension.
        png, pdf, ps, eps, and svg supported.
    frame_list : set, default = {}
        Set of frames to use in csv, if empty use all frames.
    **kwargs : dict, optional
        Additional keyword arguments to customize the plot:

        - fontsize : int, default = 10
            Font size for all text elements.
        - fig_width : float, default = 5
            Width of the figure in inches.
        - fig_height : float, default = 5
            Height of the figure in inches.
        - xaxis_fontsize : int, default = fontsize
            Font size for the theta-axis labels.
        - yaxis_fontsize : int, default = fontsize
            Font size for the rho-axis labels.

    Returns
    -------
    None

    See Also
    --------
    write_bottaro_to_csv : Creates CSV file that is inputted here

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
    >>> st.write_bottaro_to_csv(pdb_filename, output_csv_name, perspective_residue_num=perspective_residue, viewed_residue_num=viewed_residue)
    Output values written to testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv
    >>> st.visualize_two_residue_movement_scatterplot('testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv', 
    ...                                                 plot_outfile='testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_10frames_scatter.png')

    """
    bottaro_values = pd.read_csv(csv, sep=',')

    if frame_list:
        bottaro_values = bottaro_values[bottaro_values['frame'].isin(frame_list)]

    theta_values = bottaro_values['theta']
    theta_values_rad = np.radians(theta_values)

    rho_values = bottaro_values['rho_dist']

    ax = set_polar_grid(**kwargs)
    ax.scatter(theta_values_rad, rho_values, color = 'purple', s=1, alpha = 0.5)

    # Draw nucleotide ring in polar plot
    r_for_ring = np.ones(7)*1.3
    theta_for_ring = np.linspace(0, 2 * np.pi, 7)   
    ax.fill(theta_for_ring,r_for_ring, color = 'black', fill=False)

    r_for_purine = np.array([1.3, 1.3, 2.58575693, 3.075, 2.58575693, 3.49, 2.58575693, 1.3])
    theta_for_purine = np.array([4*np.pi/3, np.pi, -3.04534743, -2.61795,-2.19064032, -1.88, -2.19064032, 4*np.pi/3])  
    ax.fill(theta_for_purine, r_for_purine, color = 'black', fill=False, alpha = 0.5)

    if plot_outfile:
        plt.savefig(plot_outfile)
    else:
        plt.show()

def visualize_two_residue_movement_heatmap(csv: str, plot_outfile: str = '', frame_list: set = {}, **kwargs) -> None:
    """
    Creates heatmap of two-residue movement relative to each other.

    2D shaded contour plot of the density of points in the 
    visualize_two_residue_movement_scatterplot() scatterplot.

    Parameters
    ----------
    csv : str
        Filepath to csv file containing data on the movement
        of two residues relative to each other (r, rho, and theta values). Created
        in residue_movement.
    plot_outfile : str
        Filepath of the image file to write to. Format inferred from file extension.
        png, pdf, ps, eps, and svg supported.
    frame_list : set, default = {}
        Set of frames to use in csv, if empty use all frames.
    **kwargs : dict, optional
        Additional keyword arguments to customize the plot:
        
        - fontsize : int, default = 10
            Font size for all text elements.
        - fig_width : float, default = 5
            Width of the figure in inches.
        - fig_height : float, default = 5
            Height of the figure in inches.
        - xaxis_fontsize : int, default = fontsize
            Font size for the theta-axis labels.
        - yaxis_fontsize : int, default = fontsize
            Font size for the rho-axis labels.
        - cb_fontsize : int, default = 10
            Font size for the colorbar tick labels.

    Returns
    -------
    None

    See Also
    --------
    write_bottaro_to_csv : Creates CSV file that is inputted here
    
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
    >>> st.visualize_two_residue_movement_heatmap('testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv', 
    ...                                                 plot_outfile='testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_10frames_heat.png')

    """
    bottaro_values = pd.read_csv(csv, sep=',')

    if frame_list:
        bottaro_values = bottaro_values[bottaro_values['frame'].isin(frame_list)]

    theta_values = bottaro_values['theta']
    theta_values_rad = np.radians(theta_values)

    # convert rho values from nm to Angstroms
    rho_values = bottaro_values['rho_dist']

    ax = set_polar_grid(**kwargs)
    ax = kdeplot(x=theta_values_rad, y=rho_values, fill=True, levels = [0.1*i for i in range(1,11)], cbar = True, cmap = 'gist_earth_r')
    plt.xlabel('')
    plt.ylabel('')

    cbar = ax.collections[-1].colorbar
    levels = [0.1*i for i in range(1,11)]
    formatted_labels = ['{:.1f}'.format(level) for level in levels]
    cbar.set_ticklabels(formatted_labels)
    cbar.ax.tick_params(labelsize=kwargs.get('cb_fontsize', 10))
    cbar.ax.set_position([0.85, 0.15, 0.05, 0.7])  
    cbar.ax.set_title('Density', fontsize=kwargs.get('legend_fontsize', 10), pad=15)

    # Draw nucleotide ring in polar plot
    r_for_ring = np.ones(7)*1.3
    theta_for_ring = np.linspace(0, 2 * np.pi, 7)   
    ax.fill(theta_for_ring,r_for_ring, color = 'black', fill=False)

    r_for_purine = np.array([1.3, 1.3, 2.58575693, 3.075, 2.58575693, 3.49, 2.58575693, 1.3])
    theta_for_purine = np.array([4*np.pi/3, np.pi, -3.04534743, -2.61795,-2.19064032, -1.88, -2.19064032, 4*np.pi/3])  
    ax.fill(theta_for_purine, r_for_purine, color = 'black', fill=False, alpha = 0.5)

    if plot_outfile:
        plt.savefig(plot_outfile)
    else:
        plt.show()

@functools.wraps(visualize_two_residue_movement_scatterplot)
def display_psf_scatter(*args, **kwargs):
    return visualize_two_residue_movement_scatterplot(*args, **kwargs)

display_psf_scatter.__doc__ = f"""
Alias for `visualize_two_residue_movement_scatterplot()`.

{visualize_two_residue_movement_scatterplot.__doc__}
"""

@functools.wraps(visualize_two_residue_movement_heatmap)
def display_psf_heatmap(*args, **kwargs):
    return visualize_two_residue_movement_heatmap(*args, **kwargs)

display_psf_heatmap.__doc__ = f"""
Alias for `visualize_two_residue_movement_heatmap()`.

{visualize_two_residue_movement_heatmap.__doc__}
"""

# -----------------------------------------------------------------------------
# PSF v2 plotting helpers (amino-acid capable)
# -----------------------------------------------------------------------------

def _convex_hull_2d(points_xy: np.ndarray) -> np.ndarray:
    """Compute the convex hull of 2D points (Monotonic chain).

    Parameters
    ----------
    points_xy : (N,2) ndarray

    Returns
    -------
    hull : (M,2) ndarray
        Hull vertices in counterclockwise order.
    """
    pts = np.asarray(points_xy, dtype=float)
    if pts.shape[0] <= 2:
        return pts

    # Sort by x then y
    pts = pts[np.lexsort((pts[:, 1], pts[:, 0]))]

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(tuple(p))

    upper = []
    for p in pts[::-1]:
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(tuple(p))

    hull = lower[:-1] + upper[:-1]
    return np.array(hull, dtype=float)

def _find_atom_index_by_name(top: "md.Topology", residue_index: int, atom_name: str) -> int:
    """Return atom.index for atom_name within residue_index (0-based mdtraj index)."""
    for a in top.residue(residue_index).atoms:
        if a.name == atom_name:
            return int(a.index)
    raise ValueError(f"Could not find atom '{atom_name}' in residue index {residue_index}.")

'''
def _draw_perspective_residue_outline(ax: mpl.projections.polar.PolarAxes,
                                     pdb: str,
                                     residue_index: int = 0,
                                     frame: int = 0,
                                     ring_atoms: set | None = None,
                                     anchor_A0: str | None = None,
                                     anchor_A1: str | None = None,
                                     **kwargs) -> None:
    """Draw the perspective residue outline + anchor marker in the polar PSF plot.

    This projects the residue's ring atoms into the same in-plane basis used for
    Bottaro (rho, theta) and draws the convex hull as the outline.

    Notes
    -----
    Requires mdtraj and the residue_movement SVD helpers.
    """
    if md is None:
        raise ImportError("mdtraj is required for PSF v2 plotting (install mdtraj).")

    from . import residue_movement as rm  # local import to avoid circular imports

    trj = md.load(pdb)
    top = trj.topology

    if residue_index < 0 or residue_index >= top.n_residues:
        raise ValueError(f"residue_index={residue_index} out of range for PDB with {top.n_residues} residues")

    if frame < 0 or frame >= trj.n_frames:
        raise ValueError(f"frame={frame} out of range for PDB with {trj.n_frames} frames")

    resname = top.residue(residue_index).name

    if ring_atoms is None:
        ring_atoms = rm.determine_default_ring_atoms(resname)
    if not ring_atoms:
        return  # nothing to draw

    # Anchor defaults
    if anchor_A0 is None or anchor_A1 is None:
        if hasattr(rm, 'determine_default_anchors_v2'):
            anchor_A0, anchor_A1 = rm.determine_default_anchors(resname)
        else:
            anchor_A0, anchor_A1 = rm.determine_default_anchors(resname)

    # Collect atom indices
    ring_idx = []
    for a in top.residue(residue_index).atoms:
        if a.name in ring_atoms:
            ring_idx.append(int(a.index))
    ring_idx = np.array(ring_idx, dtype=int)
    if ring_idx.size < 3:
        return

    A0_idx = _find_atom_index_by_name(top, residue_index, anchor_A0)
    A1_idx = _find_atom_index_by_name(top, residue_index, anchor_A1)

    xyzA = trj.xyz[frame] * 10.0  # nm -> Å
    ringA = xyzA[ring_idx]
    A0A = xyzA[A0_idx]
    A1A = xyzA[A1_idx]

    # Plane + basis (use the same SVD + anchor sign fixing as PSF v2)
    C, n = rm._svd_plane_from_coords(ringA)
    n = rm._fix_normal_sign(n, C, A0A, A1A)
    xhat, yhat = rm._in_plane_basis(C, n, A0A, ringA)

    # Project ring atoms into plane basis
    pts = []
    for p in ringA:
        v = p - C
        v_par = v - np.dot(v, n) * n
        pts.append([float(np.dot(v_par, xhat)), float(np.dot(v_par, yhat))])
    pts = np.asarray(pts, dtype=float)

    hull = _convex_hull_2d(pts)
    if hull.shape[0] < 3:
        return

    # Convert hull to polar (theta, rho)
    x = hull[:, 0]
    y = hull[:, 1]
    theta = (np.arctan2(y, x) % (2 * np.pi))
    rho = np.sqrt(x * x + y * y)

    # Close the polygon
    theta = np.append(theta, theta[0])
    rho = np.append(rho, rho[0])

    ax.plot(theta, rho, color='black', linewidth=kwargs.get('outline_lw', 0.9), alpha=kwargs.get('outline_alpha', 0.9))

    # Anchor marker (A0)
    v0 = A0A - C
    v0_par = v0 - np.dot(v0, n) * n
    x0 = float(np.dot(v0_par, xhat))
    y0 = float(np.dot(v0_par, yhat))
    th0 = float(np.arctan2(y0, x0) % (2 * np.pi))
    rh0 = float(np.sqrt(x0 * x0 + y0 * y0))
    ax.scatter([th0], [rh0], s=kwargs.get('anchor_size', 18), color='black', alpha=kwargs.get('anchor_alpha', 0.9))
    ax.text(th0, rh0 + kwargs.get('anchor_text_pad', 0.25), anchor_A0,
            ha='center', va='center', fontsize=kwargs.get('anchor_fontsize', kwargs.get('fontsize', 10)))

    # Optional residue label
    if kwargs.get('show_resname', True):
        ax.text(0, 0, resname, ha='center', va='center', fontsize=kwargs.get('resname_fontsize', kwargs.get('fontsize', 10)))
'''

def _draw_perspective_residue_outline(ax,
                                      pdb: str,
                                      residue_index: int = 0,
                                      frame: int = 0,
                                      outline_lw: float = 1.2,
                                      outline_alpha: float = 0.9,
                                      anchor_marker_size: float = 18,
                                      anchor_label: bool = True,
                                      **kwargs) -> None:
    """
    Draw the *actual* perspective residue ring(s) projected into the Bottaro plane.

    - Purines (A/G/INO): draws 5-member + 6-member fused rings
    - Pyrimidines (C/U/T): draws 6-member ring
    - PHE/TYR: draws phenyl ring
    - HIS: draws imidazole
    - TRP: draws fused indole rings
    - ARG: draws guanidinium as a 4-atom diamond-ish outline (NE-CZ-NH1/NH2)
    """
    import mdtraj as md
    import numpy as np
    import math

    # These helpers should exist in residue_movement.py in your v2 work:
    # determine_default_ring_atoms(res_name) -> Set[str]
    # determine_default_anchors(res_name) -> (A0_name, A1_name)
    # _svd_plane_from_coords(coordsA) -> (C, n)
    # _fix_normal_sign(n, C, A0, A1) -> n_signed
    # _in_plane_basis(C, n, A0, ring_coordsA) -> (xhat, yhat)
    from .residue_movement import (
        determine_default_ring_atoms,
        determine_default_anchors,
        _svd_plane_from_coords,
        _fix_normal_sign,
        _in_plane_basis
    )

    trj = md.load(pdb)
    top = trj.topology
    if trj.n_frames == 0:
        return
    if frame < 0 or frame >= trj.n_frames:
        frame = 0

    residues = list(top.residues)
    if residue_index < 0 or residue_index >= len(residues):
        raise ValueError(f"residue_index={residue_index} out of range for PDB with {len(residues)} residues.")

    res = residues[residue_index]
    res_name = res.name

    ring_atomnames = determine_default_ring_atoms(res_name)
    if not ring_atomnames:
        return

    A0_name, A1_name = determine_default_anchors(res_name)

    xyzA = trj.xyz[frame] * 10.0  # nm -> Å

    # Grab ring atom coords + anchor coords
    ring_atoms = [a for a in res.atoms if a.name in ring_atomnames]
    if len(ring_atoms) < 3:
        return

    def _find_atom_in_res(atomnm: str):
        for a in res.atoms:
            if a.name == atomnm:
                return a
        raise ValueError(f"Anchor atom '{atomnm}' not found in residue {res_name} (index={residue_index}).")

    A0_atom = _find_atom_in_res(A0_name)
    A1_atom = _find_atom_in_res(A1_name)

    ring_coordsA = np.array([xyzA[a.index] for a in ring_atoms], dtype=float)
    A0A = xyzA[A0_atom.index].astype(float)
    A1A = xyzA[A1_atom.index].astype(float)

    # Plane + basis (same logic as PSF computation)
    C, n = _svd_plane_from_coords(ring_coordsA)
    n = _fix_normal_sign(n, C, A0A, A1A)
    xhat, yhat = _in_plane_basis(C, n, A0A, ring_coordsA)

    # Project atom positions into plane coords (x,y)
    pos2d = {}
    for a in ring_atoms + [A0_atom]:
        v = xyzA[a.index] - C
        x = float(np.dot(v, xhat))
        y = float(np.dot(v, yhat))
        pos2d[a.name] = (x, y)

    atomset = set(pos2d.keys())

    def _paths_for_res(resname: str, atomnames: set[str]):
        # Purines (A/G/INO): fused rings
        if {"N9", "C8", "N7", "C5", "C4"}.issubset(atomnames):
            return [
                ["N9", "C8", "N7", "C5", "C4"],              # 5-member
                ["C4", "N3", "C2", "N1", "C6", "C5"],        # 6-member
            ]
        # Pyrimidines
        if {"N1", "C2", "N3", "C4", "C5", "C6"}.issubset(atomnames):
            return [["N1", "C2", "N3", "C4", "C5", "C6"]]

        # PHE/TYR ring
        if {"CG", "CD1", "CE1", "CZ", "CE2", "CD2"}.issubset(atomnames):
            return [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]]

        # HIS ring
        if {"CG", "ND1", "CE1", "NE2", "CD2"}.issubset(atomnames):
            return [["CG", "ND1", "CE1", "NE2", "CD2"]]

        # TRP (indole): fused 5 + 6
        if {"CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"}.issubset(atomnames):
            return [
                ["CG", "CD1", "NE1", "CE2", "CD2"],                 # 5-member
                ["CD2", "CE2", "CE3", "CZ3", "CH2", "CZ2"],         # 6-member
            ]

        # ARG guanidinium: 4 atoms
        if {"NE", "CZ", "NH1", "NH2"}.issubset(atomnames):
            return [["NE", "CZ", "NH1", "NH2"]]

        # Fallback: draw nothing (or you could implement convex hull)
        return []

    paths = _paths_for_res(res_name, atomset)

    # Draw each path as a closed polygon in polar coords
    for path in paths:
        if any(a not in pos2d for a in path):
            continue
        xy = np.array([pos2d[a] for a in path], dtype=float)

        rho = np.sqrt(xy[:, 0] ** 2 + xy[:, 1] ** 2)
        theta = np.arctan2(xy[:, 1], xy[:, 0])

        # close
        rho = np.append(rho, rho[0])
        theta = np.append(theta, theta[0])

        ax.plot(theta, rho, color="black", linewidth=outline_lw, alpha=outline_alpha)

    # Anchor marker (A0)
    if A0_name in pos2d:
        x0, y0 = pos2d[A0_name]
        th0 = math.atan2(y0, x0)
        r0 = math.sqrt(x0 * x0 + y0 * y0)
        ax.scatter([th0], [r0], s=anchor_marker_size, color="black", zorder=10)
        if anchor_label:
            ax.text(th0, r0, A0_name, fontsize=kwargs.get("ring_label_fontsize", 9))

def set_polar_grid_v2(**kwargs) -> mpl.projections.polar.PolarAxes:
    """Set up axes for PSF v2 (generic, residue-agnostic)."""
    fontsize = kwargs.get('fontsize', 10)
    plt.rcParams.update({'font.size': fontsize})

    fig = plt.figure(figsize=(kwargs.get('fig_width', 5), kwargs.get('fig_height', 5)))
    ax = fig.add_subplot(111, polar=True)

    # Theta ticks
    tick_degs = kwargs.get('theta_tick_degrees', None)
    if tick_degs is None:
        tick_degs = [0, 90, 180, 270]
    tick_rads = np.radians(tick_degs)
    ax.set_xticks(tick_rads)
    ax.set_xticklabels([fr"$\theta={d}^\circ$" for d in tick_degs], fontsize=kwargs.get('xaxis_fontsize', fontsize))
    ax.tick_params(pad=-10)

    rho_max = kwargs.get('rho_max', 10)
    ax.set_rlim(0, rho_max)
    ax.set_rticks(np.linspace(0, rho_max, 3, endpoint=True))
    ax.tick_params(axis='y', labelsize=kwargs.get('yaxis_fontsize', fontsize))
    ax.set_rlabel_position(180)
    plt.text(x=np.radians(178), y=rho_max + 2, s=r"$\rho\text{ }(\AA)$",
             ha="center", va='center', fontsize=kwargs.get('fontsize', fontsize))

    ax.grid(color='gray', linestyle='--', linewidth=0.5)
    return ax

def visualize_two_residue_movement_scatterplot_v2(csv: str,
                                                 pdb: str,
                                                 plot_outfile: str = '',
                                                 frame_list: set = {},
                                                 perspective_residue_index: int = 0,
                                                 outline_frame: int = 0,
                                                 **kwargs) -> None:
    """PSF scatterplot (v2): amino-acid compatible ring overlay.

    This is a drop-in alternative to `visualize_two_residue_movement_scatterplot`.
    It uses the PSF CSV for points and the 2-residue PDB (from `filter_traj_to_pdb`)
    to draw the *actual* perspective residue outline and anchor marker.

    Parameters
    ----------
    csv : str
        CSV produced by `write_bottaro_to_csv_svd` (or legacy writer).
    pdb : str
        2-residue PDB used to create the CSV (must contain the perspective residue).
    perspective_residue_index : int, default=0
        mdtraj residue index in the PDB (0-based). For a 2-residue PDB this is usually 0.
    outline_frame : int, default=0
        Frame from the PDB to use for drawing the ring outline.

    Notes
    -----
    If you want legacy nucleotide tick labels (C2/C4/C6), use the old function.
    """

    bottaro_values = pd.read_csv(csv, sep=',')
    if frame_list:
        bottaro_values = bottaro_values[bottaro_values['frame'].isin(frame_list)]

    theta_values_rad = np.radians(bottaro_values['theta'].to_numpy())
    rho_values = bottaro_values['rho_dist'].to_numpy()

    ax = set_polar_grid_v2(**kwargs)
    ax.scatter(theta_values_rad, rho_values,
               color=kwargs.get('point_color', 'purple'),
               s=kwargs.get('point_size', 1),
               alpha=kwargs.get('point_alpha', 0.5))

    _draw_perspective_residue_outline(
        ax,
        pdb=pdb,
        residue_index=perspective_residue_index,
        frame=outline_frame,
        **kwargs
    )

    if plot_outfile:
        plt.savefig(plot_outfile)
    else:
        plt.show()

def visualize_two_residue_movement_heatmap_v2(csv: str,
                                             pdb: str,
                                             plot_outfile: str = '',
                                             frame_list: set = {},
                                             perspective_residue_index: int = 0,
                                             outline_frame: int = 0,
                                             **kwargs) -> None:
    """PSF heatmap (v2): amino-acid compatible ring overlay."""

    bottaro_values = pd.read_csv(csv, sep=',')
    if frame_list:
        bottaro_values = bottaro_values[bottaro_values['frame'].isin(frame_list)]

    theta_values_rad = np.radians(bottaro_values['theta'].to_numpy())
    rho_values = bottaro_values['rho_dist'].to_numpy()

    ax = set_polar_grid_v2(**kwargs)
    ax = kdeplot(
        x=theta_values_rad,
        y=rho_values,
        fill=True,
        levels=[0.1 * i for i in range(1, 11)],
        cbar=True,
        cmap=kwargs.get('heatmap_cmap', 'gist_earth_r')
    )
    plt.xlabel('')
    plt.ylabel('')

    # Optional: format colorbar like legacy
    try:
        cbar = ax.collections[-1].colorbar
        levels = [0.1 * i for i in range(1, 11)]
        formatted_labels = ['{:.1f}'.format(level) for level in levels]
        cbar.set_ticklabels(formatted_labels)
        cbar.ax.tick_params(labelsize=kwargs.get('cb_fontsize', 10))
        cbar.ax.set_title('Density', fontsize=kwargs.get('legend_fontsize', 10), pad=15)
        
    # TODO Remove later 
        # cbar = ax.collections[-1].colorbar
        # cbar.set_ticks([])  # remove tick marks and tick labels
        # cbar.ax.set_title('Relative\ndensity',
        #                 fontsize=kwargs.get('legend_fontsize', 10),
        #                 pad=15)
    except Exception:
        pass

    _draw_perspective_residue_outline(
        ax,
        pdb=pdb,
        residue_index=perspective_residue_index,
        frame=outline_frame,
        **kwargs
    )

    if plot_outfile:
        plt.savefig(plot_outfile)
    else:
        plt.show()


@functools.wraps(visualize_two_residue_movement_scatterplot_v2)
def display_psf_scatter_v2(*args, **kwargs):
    return visualize_two_residue_movement_scatterplot_v2(*args, **kwargs)

display_psf_scatter_v2.__doc__ = f"""
Alias for `visualize_two_residue_movement_scatterplot_v2()`.

{visualize_two_residue_movement_scatterplot_v2.__doc__}
"""


@functools.wraps(visualize_two_residue_movement_heatmap_v2)
def display_psf_heatmap_v2(*args, **kwargs):
    return visualize_two_residue_movement_heatmap_v2(*args, **kwargs)

display_psf_heatmap_v2.__doc__ = f"""
Alias for `visualize_two_residue_movement_heatmap_v2()`.

{visualize_two_residue_movement_heatmap_v2.__doc__}
"""


class NoResidues(Exception):
    """Raised if user tries to make SSF with a trajectory of <1 residue"""
    pass

if __name__ == '__main__':
    # 10 frame test
    visualize_two_residue_movement_scatterplot('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv')

    # Multi-frame test
    visualize_two_residue_movement_scatterplot('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv')

    # Multi-frame heatmap test
    visualize_two_residue_movement_heatmap('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv')

    # Write to outfile tests
    output = "stacker/testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_3200frames_scatter.png"
    create_parent_directories(output)
    visualize_two_residue_movement_scatterplot('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv', plot_outfile='testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_3200frames_scatter.png')
    visualize_two_residue_movement_heatmap('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv', plot_outfile='testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_3200frames_heatmap.png')

