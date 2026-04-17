"""
Run K Means Clustering and Principal Component Analysis

This module contains functions to run K Means Clustering on SSF
results and visualize the clusters with barplots, silhouette analysis,
and PCA scatterplots. 

"""

import os
import numpy as np
from numpy import typing
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.decomposition import PCA
import seaborn as sns


### VARIABLES ###

N_RESIDUES = 127
RANGE_N_CLUSTERS = [2,3,4,5,6,7,8] # Range of n_clusters to try for kmeans
dataset_names = ['../testing/5JUP_N2_tGGG_aCCU_+1GCU_data.txt.gz', '../testing/5JUP_N2_tGGG_aCCU_+1CGU_data.txt.gz']  # Add more dataset names as needed
filepaths = ['5JUP_N2_tGGG_aCCU_+1GCU', '5JUP_N2_tGGG_aCCU_+1CGU']  # Add more dataset names as needed
indir = '~/Downloads/stacker/testing/' # Directory with data.txt output from StACKER (created with -d flag)
outdir = '~/Downloads/stacker/testing/' # Outdir for clustering results and kmeans plot

##################

def _flatten_ssf_frames(data: np.ndarray) -> np.ndarray:
    """
    Convert SSF data to shape (n_frames, n_features).

    Accepts:
      - 3D array: (n_frames, n_res, n_res)
      - 2D array: already flattened (n_frames, n_features)
    """
    data = np.asarray(data)

    if data.ndim == 3:
        n_frames = data.shape[0]
        return data.reshape(n_frames, -1)

    if data.ndim == 2:
        return data

    raise ValueError(
        f"Expected SSF array with 2 or 3 dimensions, got shape {data.shape}"
    )


# def read_and_preprocess_data(dataset_names) -> dict:
#     """
#     read_and_preprocess_data(
#         (file1, file2, ...)
#     )

#     Reads and preprocesses SSF data for K Means analysis per dataset.

#     Reads SSF data from txt files for each dataset, decompresses the data, 
#     and attaches each Trajectory to its frame-wise SSF results. The values
#     are flattened SSF lists, so rather than a 3200 frames x 127 res x 127 res, 
#     it's a 3200 frames x 16129 res-res pairs. For example, a 2-residue, 2-frame
#     SSF of 

#         [ [[1, 2],
#         [3, 4]],
        
#         [[5, 6],
#         [7, 8]] ]

#     is flattened to:

#         [[1, 2, 3, 4],
#         [5, 6, 7, 8]]

#     Parameters
#     ----------
#     file1, file2, ... : list of str
#         List of filenams to read and preprocess.
#         Outputted from `-s ssf -d output.txt`.
#         Should be in the format {datapath}/{traj_name}.txt.gz
        
#     Returns
#     -------
#     data_arrays : dict
#         Dictionary where keys are dataset names and values are the processed data arrays.

#     See Also
#     --------
#     create_kmeans_input : Stacks SSF data into a single 2D Numpy array.    
    
#     Examples
#     --------
#     >>> import stacker as st
#     >>> dataset_names = ['testing/5JUP_N2_tGGG_aCCU_+1GCU.txt.gz', 'testing/5JUP_N2_tGGG_aCCU_+1CGU.txt.gz']  # 3200 frames, SSFs of 127 x 127 residues
#     >>> data_arrays = st.read_and_preprocess_data(dataset_names)
#     >>> print(data_arrays['dataset1'].shape)
#     (3200, 16129)

#     """
#     data_arrays = {}
#     for filepath in dataset_names:
#         file = filepath.split('/')[-1]
#         name = file.split('.')[0]

#         print('Reading data:', file)
#         data = np.loadtxt(filepath)
#         data_arrays[name] = data
#     return data_arrays

def read_and_preprocess_data(dataset_names: list[str]) -> dict[str, np.ndarray]:
    """
    Read SSF datasets and return flattened frame-wise arrays.

    Supports:
      - .npy / .npz via np.load
      - .txt / .txt.gz via np.loadtxt
    """
    data_arrays = {}

    for filepath in dataset_names:
        file = os.path.basename(filepath)
        name = file.split('.')[0]
        print("Reading data:", file)

        if filepath.endswith(".npy"):
            raw = np.load(filepath)

        elif filepath.endswith(".npz"):
            loaded = np.load(filepath)
            if "arr_0" in loaded:
                raw = loaded["arr_0"]
            else:
                keys = list(loaded.keys())
                if len(keys) != 1:
                    raise ValueError(
                        f"{filepath} has multiple arrays; expected one. Keys: {keys}"
                    )
                raw = loaded[keys[0]]

        else:
            raw = np.loadtxt(filepath)

        data_arrays[name] = _flatten_ssf_frames(raw).astype(np.float32, copy=False)

    return data_arrays

def create_kmeans_input(data_arrays: dict) -> typing.ArrayLike:
    """
    Blinds SSF Data (removes trajectory labels) for input to K Means

    Stacks SSF data into a single 2D numpy array from all frames of 
    all trajectories without labels for each frame. Used for input to
    KMeans Clustering

    Parameters
    ----------
    data_arrays : dict
        Output of read_and_preprocess_data(). Dictionary where keys are dataset 
        names and values are the processed data arrays.

    Returns
    -------
    blinded_data : np.typing.ArrayLike
        A 2D numpy array containing all frames stacked together.

    See Also
    --------
    read_and_preprocess_data : Reads and preprocesses data for each dataset

    Examples
    --------
    >>> import stacker as st
    >>> data_arrays = {
    ...     'dataset1': np.random.rand(3200, 16129),
    ...     'dataset2': np.random.rand(3200, 16129)
    ... }
    >>> kmeans_input = st.create_kmeans_input(data_arrays)
    >>> print(kmeans_input.shape)
    (6400, 16129)
    """
    blind_data = list(data_arrays.values())
    data = np.vstack(blind_data)
    print(data.shape)
    return data

def run_kmeans(data_arrays : dict, n_clusters: int,
               max_iter: int = 1000, n_init: int = 20, random_state: int = 1, outdir: str = '') -> np.ndarray :
    """
    Performs KMeans clustering on blinded SSF data and saves the results.

    This function applies the KMeans clustering algorithm to the provided
    blinded SSF data, assigns each frame to a cluster, and counts the number of
    frames in each cluster for each dataset. The results are printed and
    saved to a file.

    Parameters
    ----------
    data_arrays : dict
        Output of read_and_preprocess_data(). Dictionary where keys are dataset 
        names and values are the processed data arrays.
    n_clusters : int
        The number of clusters to form 
    max_iter : int, default=1000
        Maximum number of iterations of the k-means algorithm for a single run.
    n_init : int, default=20
        Number of times the k-means algorithm will be run with different centroid seeds.
    random_state : int, default=1
        Determines random number generation for centroid initialization.
    outdir : str, default=''
        Directory to save the clustering results.
        If empty, just prints to standard output.

    Returns
    -------
    np.ndarray
        The labels of the clusters for each frame.

    See Also
    --------
    create_kmeans_input : blinds SSF Data for input to K Means
    read_and_preprocess_data : reads and preprocesses SSF data for K Means analysis per dataset

    Examples
    --------
    >>> import stacker as st
    >>> data_arrays = {
    ...     'dataset1': np.random.rand(3200, 16129),
    ...     'dataset2': np.random.rand(3200, 16129)
    ... }
    >>> blinded_data = st.create_kmeans_input(data_arrays)
    >>> st.run_kmeans(blinded_data, n_clusters=4)
    Reading data: dataset1
    Reading data: dataset2
    (6400, 16129)
    {'dataset1': array([800, 800, 800, 800]), 'dataset2': array([800, 800, 800, 800])}
    Dataset: dataset1
        Cluster 1: 800 matrices
        Cluster 2: 800 matrices
        Cluster 3: 800 matrices
        Cluster 4: 800 matrices
    Dataset: dataset2
        Cluster 1: 800 matrices
        Cluster 2: 800 matrices
        Cluster 3: 800 matrices
        Cluster 4: 800 matrices

    """
    global blindframes_labelled_by_cluster
    global silhouette_avg
    global sample_silhouette_values

    if outdir and not outdir.endswith('/'):
        outdir += '/'

    blinded_data = create_kmeans_input(data_arrays)

    kmeans_func_instance = KMeans(n_clusters=n_clusters, max_iter=max_iter, n_init=n_init, random_state=random_state)
    blindframes_labelled_by_cluster = kmeans_func_instance.fit_predict(blinded_data)
    silhouette_avg = silhouette_score(blinded_data, blindframes_labelled_by_cluster)

    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )

    sample_silhouette_values = silhouette_samples(blinded_data, blindframes_labelled_by_cluster)

    counts = {}
    labels = blindframes_labelled_by_cluster
    for name, arr in data_arrays.items():
        counts[name] = np.bincount(labels[:len(arr)], minlength=n_clusters)
        labels = labels[len(arr):]  # Move to the next dataset

    # Print the results
    for name, count in counts.items():
        print(f'Dataset: {name}')
        for cluster in range(n_clusters):
            print(f'\tCluster {cluster+1}: {count[cluster]} matrices')

    # Save results to file
    if outdir:
        outfile_path = outdir + 'clustering_results_' + str(n_clusters) + '.txt'
        outfile = open(outfile_path, 'w')
        outfile.write('cluster trj number\n')
        for name, count in counts.items():
            for cluster in range(n_clusters):
                outfile.write(f'{cluster+1} {name} {count[cluster]}\n')
        outfile.close()
        print(f"Results written to: {outfile_path}")

    return blindframes_labelled_by_cluster

def plot_cluster_trj_data(cluster_file: str, outfile: str, x_labels_map: dict = None) -> None:
    """
    Plots the output of run_kmeans() to a PNG file.

    Creates a grouped bar plot of the number of frames from each trajectory in each cluster
    following KMeans clustering. Writes the plot output to a PNG file.

    Parameters
    ----------
    input_file : str
        Path to clustering results written by run_kmeans()
    outfile : str
        Filepath where the plot PNG file will be saved.
    x_labels_map : dict, optional
        Dictionary to remap x labels. Keys are original labels and values are new labels.
        
    Returns
    -------
    None

    Examples
    --------
    This will read the clustering results from 'clustering_results.txt',
    create a bar plot, and save it as 'kmeans_plot.cluster_4.png' in the 
    specified output directory.

    >>> import stacker as st
    >>> st.plot_cluster_trj_data('clustering_results.txt', "../testing/kmeans_plot.png", {'5JUP_N2_tGGG_aCCU_+1CGU_data': 'tGGG_aCCU_+1CGU', '5JUP_N2_tGGG_aCCU_+1GCU_data': 'tGGG_aCCU_+1GCU'})

    """
    cluster_data = pd.read_table(cluster_file, sep=' ', header=0, quotechar="\"")
    sns.set_theme(style="white", font_scale=1.2)

    g = sns.FacetGrid(cluster_data.dropna(subset=['trj']), col="cluster", col_wrap=2, height=6, despine=False)
    colors = sns.color_palette("deep", len(cluster_data['trj'].unique()))
    g.map(plt.bar, 'trj', 'number', color=colors) 

    for ax in g.axes.flat:
        if x_labels_map:
            # Get unique trajectory names
            unique_trjs = cluster_data['trj'].unique()
            labels = [x_labels_map.get(trj, trj) for trj in unique_trjs]
            print(f"Original labels: {unique_trjs}")
            print(f"Mapped labels: {labels}")
            ax.set_xticks(range(len(labels)))  # Set fixed ticks
            ax.set_xticklabels(labels)
        
        for label in ax.get_xticklabels():
            label.set_rotation(90)
            label.set_ha('right')

        ax.set_xlabel("Trajectory")
        ax.set_ylabel("Number of Frames")  

    for i, ax in enumerate(g.axes.flat):
        ax.set_title(f"Cluster {i + 1}")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f"Plot Outputted to {outfile}")
    plt.close()

def plot_silhouette(n_clusters : int, blind_data : typing.ArrayLike, outdir : str = ''):
    '''
    Creates Silhouette plots to determine the best number of clusters

    Parameters
    ----------
    n_clusters : int, default = 0
        The number of clusters to form.
    blind_data : np.typing.ArrayLike
        A 2D numpy array containing all frames stacked together.
        Output of create_kmeans_input()
    outfile : str
        Filepath where the plot PNG file will be saved.
    '''
    if outdir and not outdir.endswith('/'):
        outdir += '/'

    plt.figure(figsize=(10, 7))
    plt.xlim([-1, 1])
    plt.ylim([0, len(blind_data) + (n_clusters + 1) * 10])

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[blindframes_labelled_by_cluster == i]
        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        plt.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        plt.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    plt.title("The silhouette plot for the various clusters.")
    plt.xlabel("The silhouette coefficient values")
    plt.ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    plt.axvline(x=silhouette_avg, color="red", linestyle="--")

    plt.yticks([])  # Clear the yaxis labels / ticks
    plt.xticks(np.arange(-1, 1.1, 0.1))

    plt.tight_layout()
    plot_outpath = f"{outdir}silhouette{n_clusters}.png"
    plt.savefig(plot_outpath)
    print(f"File saved to: {plot_outpath}")
    plt.close()

"""
def plot_pca(blinded_data: np.ndarray, dataset_names: list,
             coloring: str = 'dataset', outdir: str = '', 
             cluster_labels: np.ndarray = None, new_dataset_names: dict = None) -> None:
    '''
    Creates PCA Plot to compare systems in 2D 

    Creates a PCA plot that can be colored by the KMeans clustering result
    or by dataset. Compares SSFs similarly to K Means.

    Parameters
    ----------
    blinded_data : np.ndarray
        A 2D numpy array containing all frames stacked together.
        Output of create_kmeans_input()
    dataset_names : list of str
        List of filenames to read and preprocess.
        Outputted from `stacker -s ssf -d output.txt.gz`.
        Should be in the format {datapath}/{traj_name}.txt.gz
    coloring : {'dataset', 'kmeans', 'facet'}
        Method to color the points on the scatterplot. Options:
        - dataset:  Plot all points on the same scatterplot and color by dataset of origin.
        - kmeans: Plot all points on the same scatterplot and color by KMeans Cluster with n_clusters
        - facet: Same as dataset but plot each dataset on a different coordinate grid.
    outdir : str, default=''
        Directory to save the clustering results.
    cluster_labels : np.ndarray, optional
        The labels of the clusters for each frame, output from run_kmeans.
        Used if coloring = "kmeans" to color points by cluster
    new_dataset_names : dict, optional
        Dictionary to remap dataset names. Keys are original filenames in ``dataset_names`` and values are shortened names.
        
    Returns
    -------
    None

    See Also
    --------
    create_kmeans_input : blinds SSF Data for input to K Means
    read_and_preprocess_data : reads and preprocesses SSF data for K Means analysis per dataset
    sklearn.decomposition.PCA : Runs PCA
    
    '''
    if outdir and not outdir.endswith('/'):
        outdir += '/'

    if new_dataset_names:
        dataset_names = [new_dataset_names[filepath] for filepath in dataset_names]
    else:
        dataset_names = [filepath.split('/')[-1].split('.')[0] for filepath in dataset_names]

    n_datasets = len(dataset_names)
    pca = PCA(n_components=2)
    data_reduced = pca.fit_transform(blinded_data)
    colors = []
    section_size = data_reduced.shape[0] // n_datasets
    for i in range(n_datasets):
        colors.extend([dataset_names[i]] * section_size)

    df = pd.DataFrame({
        'Principal Component 1': data_reduced[:, 0],
        'Principal Component 2': data_reduced[:, 1],
        'Color': colors
    })

    if coloring == 'facet':
        g = sns.FacetGrid(df, col='Color', col_wrap=2, height=4, despine=False, hue="Color")
        g.map_dataframe(sns.scatterplot, x='Principal Component 1', y='Principal Component 2', linewidth = 0, s = 10)


        g.set_titles(col_template='{col_name}')
        g.set_axis_labels('Principal Component 1', 'Principal Component 2')


        outfile = f"{outdir}pca_plot.by_facet.png"
        plt.savefig(outfile)
        plt.close()
    elif coloring == 'dataset':
        plt.figure(figsize=(10, 7))
        unique_colors = {name: idx for idx, name in enumerate(dataset_names)}
        df['Color'] = df['Color'].map(unique_colors)
        
        cmap = plt.get_cmap('tab10', len(dataset_names))

        scatter = plt.scatter(df['Principal Component 1'], df['Principal Component 2'], c=df['Color'], cmap=cmap, s=10)
        
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(i), markersize=10) for i in range(len(dataset_names))]
        plt.legend(handles, dataset_names, title='Dataset', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.title('PCA-reduced data by dataset')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        
        outfile = f"{outdir}pca_plot.by_dataset.png"
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
    elif coloring == "kmeans" and cluster_labels is not None:
        df['Color'] = cluster_labels
        plt.figure(figsize=(10, 7))
        
        unique_clusters = np.unique(cluster_labels)
        cmap = plt.get_cmap('tab10', len(unique_clusters))
        
        scatter = plt.scatter(df['Principal Component 1'], df['Principal Component 2'], c=df['Color'], cmap=cmap, s=10)
        
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(i), markersize=10) for i in range(len(unique_clusters))]
        plt.legend(handles, unique_clusters, title='Cluster Label', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.title('PCA-reduced data with KMeans clustering')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        
        n_clusters = np.unique(cluster_labels).size
        outfile = f"{outdir}pca_plot{n_clusters}by_cluster.png"
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
    else:
        print(f"{coloring} not a supported coloring")
        return None
"""

def plot_pca(blinded_data: np.ndarray, dataset_names: list,
             coloring: str = 'dataset', outdir: str = '',
             cluster_labels: np.ndarray = None, new_dataset_names: dict = None,
             data_arrays: dict | None = None) -> None:

    if outdir and not outdir.endswith('/'):
        outdir += '/'

    if new_dataset_names:
        display_names = [new_dataset_names.get(fp, os.path.basename(fp).split('.')[0]) for fp in dataset_names]
    else:
        display_names = [os.path.basename(fp).split('.')[0] for fp in dataset_names]

    pca = PCA(n_components=2)
    data_reduced = pca.fit_transform(blinded_data)

    if data_arrays is not None:
        colors = []
        for original_fp, display_name in zip(dataset_names, display_names):
            key = os.path.basename(original_fp).split('.')[0]
            n_rows = len(data_arrays[key])
            colors.extend([display_name] * n_rows)
    else:
        # fallback: assume equal sizes
        n_datasets = len(display_names)
        section_size = data_reduced.shape[0] // n_datasets
        colors = []
        for i in range(n_datasets):
            colors.extend([display_names[i]] * section_size)

    df = pd.DataFrame({
        'Principal Component 1': data_reduced[:, 0],
        'Principal Component 2': data_reduced[:, 1],
        'Color': colors
    })

    explained = pca.explained_variance_ratio_

if __name__ == "__main__":
    data_arrays = read_and_preprocess_data(dataset_names)
    blinded_data = create_kmeans_input(data_arrays)
    plot_pca(blinded_data, 'dataset')
    plot_pca(blinded_data, 'facet')
    for N_CLUSTERS in RANGE_N_CLUSTERS:
        run_kmeans(data_arrays, N_CLUSTERS, outdir = outdir)
        cluster_file = outdir + 'clustering_results_' + str(N_CLUSTERS) + '.txt'
        outfile = f"{outdir}kmeans_plot.cluster_{N_CLUSTERS}.png"
        plot_cluster_trj_data(cluster_file, outfile = outfile)
        plot_silhouette(N_CLUSTERS, blinded_data)
        # plot_pca(blinded_data, N_CLUSTERS, 'kmeans') # Works best with only two systems