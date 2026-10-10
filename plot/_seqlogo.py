import logomaker
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple, Dict


def seqlogo(
    protein_sequences,
    proxy: str = 'sequences',   
    clustalw_path: Optional[str] = None,  
    color_scheme: str = 'chemistry',  
    show_spines: bool = False, 
    ax: Optional[plt.Axes] = None,  
    save: Optional[str] = None  
):
    """
    Generates a sequence logo from aligned protein sequences.

    Parameters
    ----------
    protein_sequences : list
        List of protein sequences to align and generate the logo from.

    proxy : str, optional
        The proxy method for the alignment. Default is 'sequences'.

    clustalw_path : Optional[str], optional
        The path to the ClustalW executable. If None, the default ClustalW path is used.

    color_scheme : str, optional
        The color scheme to use for the sequence logo. Default is 'chemistry'.

    show_spines : bool, optional
        Whether to show spines (axis lines) on the logo. Default is False.

    ax : Optional[plt.Axes], optional
        The matplotlib Axes on which to plot the logo. If None, a new figure and axes will be created.

    save : Optional[str], optional
        If provided, the plot will be saved to this path.

    Returns
    -------
    logomaker.Logo
        A Logomaker Logo object representing the sequence logo.
    """

    sequences = clustalw_aln(protein_sequences, proxy=proxy, clustalw_path=clustalw_path)
    kwag = {}
    if ax is not None:
        kwag['ax'] = ax

    logo_df = logomaker.alignment_to_matrix(sequences=sequences)
    pos_values = np.sum(logo_df,axis = 1)
    logo_df = logo_df[pos_values > max(pos_values)/3].reset_index(drop=True)
    logo = logomaker.Logo(
        logo_df,
        color_scheme=color_scheme,
        show_spines=show_spines,
        **kwag
    )

    logo.style_xticks(rotation=90, fmt="%d", anchor=0)
    if save is not None:
        logo.ax.figure.savefig(save)

    return logo

from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import shutil

def write_fasta(sequences,seuences_id = None, file_name="sequences.fasta"):
    with open(file_name, "w") as fasta_file:
        if seuences_id == None:
            for i, seq in enumerate(sequences):
                fasta_file.write(f">seq{i+1}\n{seq}\n")
        else:
            for i, seq in zip(seuences_id,sequences):
                fasta_file.write(f">{i}\n{seq}\n")
    return file_name

def clustalw_aln(input_sequences,proxy,clustalw_path,seuences_id = None,reture_alignmeant = False):
    aln_file = f"{proxy}.aln"
    fasta_file = write_fasta(input_sequences,seuences_id = seuences_id,file_name = f'{proxy}.fasta')
    if clustalw_path == None:
        clustalw_path = shutil.which("clustalw2")
    if clustalw_path:
        pass
    else:
        raise KeyError("ClustalW2 not found. Please specify the full path to clustalw2")

    clustalw_cline = ClustalwCommandline(clustalw_path, infile=fasta_file)
    stdout, stderr = clustalw_cline()

    alignment = AlignIO.read(aln_file, "clustal")
    sequences = [str(record.seq) for record in alignment]
    if reture_alignmeant:
        return alignment
    else:
        return sequences