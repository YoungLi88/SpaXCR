# spaxcr Documentation

**spaxcr** is a Python toolkit for spatial immune repertoire analysis. It links TCR/BCR clonotype results with spatial transcriptomics coordinates, stores the result as an `AnnData` object, and provides utilities for clone querying, expansion analysis, diversity analysis, spatial smoothing, density estimation, spatial visualization, and CDR3 similarity analysis.

```python
import spaxcr
```

```{toctree}
:maxdepth: 2
:caption: User Guide

tutorial
api
```

## Main Features

- Build a spatial immune repertoire `AnnData` object from MiXCR outputs or custom clone tables.
- Integrate TCR/BCR clone information with spatial coordinates, annotations, and 10x Visium image metadata.
- Query clone-to-spot, CDR3-to-spot, and spot-to-clone relationships.
- Analyze clone expansion, group intersections, repertoire diversity, spatial density, and spatial smoothing.
- Visualize clone distributions, CDR3 length distributions, spatial clone localization, KDE overlays, boundary distance rankings, CDR3 similarity heatmaps, radial trees, and sequence logos.

## Quick Start

Install the package:

```bash
pip install spaxcr
```

Import the package:

```python
import spaxcr
```

Load MiXCR clone data:

```python
adata = spaxcr.pp.load_xcr(
    clone_file="path/to/mixcr.tsv",
    cloneId_barcode_file="path/to/cloneId_barcode.tsv",
    focus_functional_clone=True,
    focus_chain=None,
    from_mixcr=True,
    add_mutation=True,
)
```

Add spatial coordinates:

```python
adata = spaxcr.pp.addloc(
    adata,
    loc_info=rna_adata.obs[["x", "y"]],
    basis="spatial",
)
```

Create a simple spatial clone plot:

```python
spaxcr.pl.spatial(
    adata,
    clones="clone_700",
    base="anno",
)
```

## Input Overview

spaxcr expects two core immune repertoire inputs:

- A clone table, either from MiXCR or a custom table containing at least `cloneId`, `topChains`, and `aaSeqCDR3`.
- A clone-barcode table with either `barcode, cloneId` or `barcode, cloneId, count`.

Spatial coordinates are added separately through a `pandas.DataFrame` with barcode index and `x`, `y` columns.

## Where to Go Next

- Start with the [Tutorial](tutorial.md) for an end-to-end workflow.
- Use the [API Parameters](api.md) page when you need detailed function arguments.
