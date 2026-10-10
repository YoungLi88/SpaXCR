# SpaXCR

[![PyPI version](https://img.shields.io/pypi/v/spaxcr.svg)](https://pypi.org/project/spaxcr/)
[![Documentation Status](https://readthedocs.org/projects/spaxcr/badge/?version=latest)](https://spaxcr.readthedocs.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

SpaXCR is a Python toolkit for spatial single-cell TCR/BCR repertoire analysis. It connects clonotype data from MiXCR or custom clone tables with spatial transcriptomics coordinates, stores the result in an `AnnData` object, and provides tools for clone querying, expansion analysis, diversity estimation, spatial smoothing, density analysis, spatial visualization, and CDR3 similarity analysis.

## Highlights

- Build a spatial immune repertoire `AnnData` object from MiXCR outputs or custom clone tables.
- Integrate TCR/BCR clones with spatial coordinates, RNA annotations, and 10x Visium image metadata.
- Query clone-to-spot, CDR3-to-spot, and spot-to-clone relationships.
- Analyze clone expansion, group intersections, repertoire diversity, spatial density, and spatial smoothing.
- Visualize clone distributions, CDR3 length distributions, spatial clone localization, KDE overlays, boundary-distance rankings, CDR3 similarity heatmaps, radial trees, and sequence logos.
- Test spatial pairing between BCR heavy/light chains or TCR beta/alpha chains.

## Documentation

The full documentation, including a step-by-step tutorial and API reference, is hosted on Read the Docs:

- [SpaXCR documentation](https://spaxcr.readthedocs.io/en/latest/)

## Installation

SpaXCR requires Python 3.8 or newer.

```bash
pip install spaxcr
```

For the latest development version:

```bash
pip install git+https://github.com/YoungLi88/SpaXCR.git
```

## Quick start

```python
import scanpy as sc
import spaxcr

clone_file = "path/to/mixcr.tsv"
clone_barcode_file = "path/to/cloneId_barcode.tsv"

adata = spaxcr.pp.load_xcr(
    clone_file=clone_file,
    cloneId_barcode_file=clone_barcode_file,
    focus_functional_clone=True,
    focus_chain=None,
    from_mixcr=True,
    add_mutation=True,
)
```

Add spatial coordinates from an RNA `AnnData` object:

```python
rna_adata = sc.read_h5ad("path/to/rna.h5ad")
rna_adata.obs[["x", "y"]] = rna_adata.obsm["spatial"]

adata = spaxcr.pp.addloc(
    adata,
    loc_info=rna_adata.obs[["x", "y"]],
    basis="spatial",
)
```

Filter clones and create a spatial plot:

```python
adata = spaxcr.pp.filter_clone(
    adata,
    min_cdr3_length=1,
    max_cdr3_length=30,
    min_cells=1,
    min_count=1,
)

spaxcr.pl.clone_distribution(adata, kind="bar")
spaxcr.pl.spatial(adata, clones="clone_700", base="anno")
```

## Main submodules

- `spaxcr.pp`: preprocessing, data loading, and clone filtering.
- `spaxcr.tl`: analysis tools such as `query`, `expansion`, `intersection`, `diversity`, `smooth`, `density`, `aggregate`, `obsrank`, `spatial_pairing`, and `XCR_pairing`.
- `spaxcr.pl`: visualization tools such as `spatial`, `expansion`, `similarity_heatmap`, `similarity_radialtree`, `seqlogo`, and distribution plots.
- `spaxcr.extra`: additional helpers, including `radialTreee`.

## Input files

SpaXCR expects two core immune repertoire inputs:

- A clone table from MiXCR, or a custom table containing at least `cloneId`, `topChains`, and `aaSeqCDR3`.
- A clone-barcode table with either `barcode, cloneId` or `barcode, cloneId, count`.

Spatial coordinates are added separately through a `pandas.DataFrame` with barcode index and `x`, `y` columns.

For detailed input requirements, see the [documentation](https://spaxcr.readthedocs.io/en/latest/).

## Dependencies

Core dependencies include:

```text
anndata
scanpy
pandas
numpy
scipy
scikit-learn
matplotlib
seaborn
tqdm
statsmodels
venn
upsetplot
logomaker
biopython
python-Levenshtein
```

`spaxcr.pl.seqlogo` additionally requires ClustalW2. It can be available on your system `PATH`, or you can pass the executable path through `clustalw_path`.

## Development

To build the documentation locally:

```bash
python -m pip install -r docs/requirements.txt
make docs
```

The generated HTML is written to `docs/_build/html/index.html`.

## License

This project is distributed under the MIT License. See [LICENSE](LICENSE) for details.
