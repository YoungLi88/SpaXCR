# spaxcr Tutorial

This tutorial walks through the main workflow used in the example notebook:

1. Load immune repertoire data.
2. Add spatial coordinates.
3. Filter clones.
4. Generate summary plots.
5. Query clones and CDR3 sequences.
6. Run expansion, diversity, smoothing, similarity, and spatial pairing analyses.

Use the paths in this tutorial as placeholders and replace them with your own files.

## Navigation

- [Project Overview](./index.md)
- [API Parameters](./api.md)

## 1. Prepare the Environment

```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import spaxcr

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.sans-serif"] = "Arial"
sc.set_figure_params(dpi=50, dpi_save=300, vector_friendly=True)
```

## 2. Load MiXCR Clone Data

MiXCR output contains clone-level sequence and chain information. The clone-barcode file connects each clone to the spatial barcode or cell barcode where it was detected.

```python
mixcr_file = "path/to/sample_mixcr.tsv"
clone_barcode_file = "path/to/sample_cloneId_barcode.tsv"

adata = spaxcr.pp.load_xcr(
    clone_file=mixcr_file,
    cloneId_barcode_file=clone_barcode_file,
    focus_functional_clone=True,
    focus_chain=None,
    from_mixcr=True,
    add_mutation=True,
)
```

Recommended settings:

- `focus_functional_clone=True`: remove non-functional CDR3 sequences containing `*` or `_`.
- `focus_chain=None`: keep all supported TCR/BCR chains, including `TRA`, `TRB`, `IGH`, `IGK`, and `IGL`.
- `focus_chain="TCR"`: keep only `TRA` and `TRB`.
- `focus_chain="BCR"`: keep only `IGH`, `IGK`, and `IGL`.
- `from_mixcr=True`: parse MiXCR fields.
- `add_mutation=True`: extract mutation information from MiXCR alignment fields.

## 3. Load a Custom Clone Table

If the clone table is not a MiXCR output, it must contain at least these columns:

```text
cloneId    topChains    aaSeqCDR3
```

Then load it with `from_mixcr=False`:

```python
clone_file = "path/to/cloneinfo.tsv"
clone_barcode_file = "path/to/clonebarcode.tsv"

adata = spaxcr.pp.load_xcr(
    clone_file=clone_file,
    cloneId_barcode_file=clone_barcode_file,
    focus_functional_clone=True,
    focus_chain=None,
    from_mixcr=False,
    add_mutation=False,
)
```

If you also want mutation information from a custom clone table, provide these additional columns:

```text
V_align    J_align    nSeqCDR3    mutation_status
```

Then set `add_mutation=True`.

## 4. Add Spatial Coordinates

The most common pattern is to reuse coordinates from an RNA `AnnData` object:

```python
rna_adata = sc.read_h5ad("path/to/rna.h5ad")
rna_adata.obs[["x", "y"]] = rna_adata.obsm["spatial"]
loc_info = rna_adata.obs[["x", "y"]]

adata = spaxcr.pp.addloc(adata, loc_info, basis="spatial")
```

If RNA barcodes contain sample suffixes but the clone barcodes do not, normalize the names before calling `addloc`:

```python
rna_adata.obs_names = rna_adata.obs_names.map(lambda x: x[:-2])
```

You can also copy RNA annotations and Visium image metadata into the XCR object:

```python
adata.obs["anno"] = rna_adata.obs["anno"]
adata.uns["anno_colors"] = rna_adata.uns.get("anno_colors")
adata.uns["spatial"] = rna_adata.uns.get("spatial")
```

## 5. Filter Clones

Filter clones by CDR3 length, count, and the number of spots/cells in which they appear:

```python
adata = spaxcr.pp.filter_clone(
    adata,
    min_cdr3_length=1,
    max_cdr3_length=30,
    min_cells=1,
    min_count=1,
)
```

Common filtering choices:

- Use `min_cdr3_length` and `max_cdr3_length` to remove suspicious CDR3 lengths.
- Use `min_cells` to keep clones detected in at least a given number of spots/cells.
- Use `min_count` to filter by a clone-level count column.
- Use `count_column` if the clone table contains multiple count-like columns.

Inspect the parsed clone table:

```python
clonedf = adata.uns["immu_repertoire"]["clone_info"]
clonedf.head()
```

Save the processed object:

```python
adata.write("sample.xcr.h5ad")
```

## 6. Basic Summary Plots

### CDR3 Length Distribution

```python
spaxcr.pl.cdr3_length_histogram(adata, chains=None)
```

Plot selected chains:

```python
spaxcr.pl.cdr3_length_histogram(
    adata,
    chains=["TRA", "IGK"],
    colors="Set2",
)
```

### V/D/J/C Gene Usage

```python
spaxcr.pl.gene_distribution(
    adata,
    region="V",
    ntop=5,
)
```

Control chain order:

```python
spaxcr.pl.gene_distribution(
    adata,
    region="J",
    ntop=4,
    order=["TRA", "TRB", "IGK"],
)
```

### Clone or Chain Distribution

Bar plot:

```python
spaxcr.pl.clone_distribution(
    adata,
    isotype=False,
    focus_cdr3=False,
    kind="bar",
)
```

Pie plot by BCR isotype:

```python
spaxcr.pl.clone_distribution(
    adata,
    isotype=True,
    kind="pie",
    fontsize=20,
    pie_width=1,
)
```

## 7. Clone or CDR3 Intersection Between Groups

Use a grouping column in `adata.obs`, such as tissue annotation, cluster label, or region:

```python
spaxcr.pl.intersection(
    adata,
    groupby="anno",
    chain=None,
    kind="venn",
)
```

Focus on one chain:

```python
spaxcr.pl.intersection(
    adata,
    groupby="anno",
    chain="TRA",
    kind="venn",
)
```

Use an UpSet plot when there are many groups:

```python
spaxcr.pl.intersection(
    adata,
    groupby="anno",
    chain="TRA",
    kind="upset",
)
```

Return a table instead of a plot:

```python
inter_df = spaxcr.tl.intersection(
    adata,
    groupby="anno",
    chain="IGH",
    focus_cdr3=False,
    mode="total",
)
```

## 8. Spatial Visualization

### Plot a Single Clone

```python
spaxcr.pl.spatial(
    adata,
    clones="clone_1",
)
```

### Plot Multiple Clones

```python
spaxcr.pl.spatial(
    adata,
    clones=["clone_1", "clone_2"],
    colors="red",
)
```

### Plot All BCR Clones

```python
bcr_clones = adata[:, adata.var["Chain"].isin(["IGL", "IGK", "IGH"])].var_names

spaxcr.pl.spatial(
    adata,
    clones=bcr_clones,
    colors="red",
    spotsize=20,
)
```

### Plot by CDR3 Sequence

```python
spaxcr.pl.spatial(
    adata,
    cdr3s="CISRDRSGDRWVF",
    colors="blue",
)
```

### Plot on a 10x Visium Image

This requires `adata.uns["spatial"]` to contain Visium image metadata.

```python
spaxcr.pl.spatial(
    adata,
    clones=bcr_clones,
    colors="red",
    spotsize=20,
    base="image",
)
```

### Plot on an Annotation Background

```python
spaxcr.pl.spatial(
    adata,
    clones=bcr_clones,
    colors="blue",
    spotsize=200,
    scatter_spotsize=1,
    base="anno",
    base_alpha=1,
)
```

### KDE Overlay

```python
spaxcr.pl.spatial(
    adata,
    clones=bcr_clones,
    colors="blue",
    base="anno",
    kde=True,
)
```

Show only KDE without clone scatter points:

```python
spaxcr.pl.spatial(
    adata,
    clones=bcr_clones,
    colors="blue",
    no_scatter=True,
    base="anno",
    kde=True,
)
```

## 9. Clone Expansion

Plot the top expanded clones for one chain:

```python
spaxcr.pl.expansion(
    adata,
    chain="IGK",
    focus_cdr3=False,
    kind="bar",
)
```

Group by CDR3 instead of cloneId:

```python
spaxcr.pl.expansion(
    adata,
    chain="IGK",
    focus_cdr3=True,
    kind="bar",
)
```

Plot expansion categories as a pie chart:

```python
spaxcr.pl.expansion(
    adata,
    chain="IGK",
    focus_cdr3=False,
    kind="pie",
    groups=(1, 10),
    labels=["Not expanded", "Expanded", "Hyperexpanded"],
)
```

Return the expansion table:

```python
expansion_df = spaxcr.tl.expansion(
    adata,
    chain="IGH",
    focus_cdr3=False,
    groups=(1, 10),
    labels=["Not expanded", "Expanded", "Hyperexpanded"],
)
```

## 10. Spatial Boundary and Distance Ranking

First detect the boundary between an inner region and an outer region:

```python
spaxcr.pl.hull_detect(
    adata,
    inner=["tumor", "tumor_nec"],
    outter=["stroma", "stroma_immu"],
    groupby="anno",
    min_distance=None,
)
```

The function writes results to:

```python
adata.uns["immu_repertoire"]["hull"]
adata.obs["distance"]
```

Then rank clones or CDR3s by distance to the boundary:

```python
spaxcr.pl.rank_distance(
    adata,
    focus_cdr3=False,
    clone_type=None,
    ntop=10,
    figsize=(10, 5),
)
```

## 11. CDR3 Similarity Analysis

Prepare a CDR3 table:

```python
cdr3sdf = adata.var[["CDR3", "Chain"]].drop_duplicates().dropna()
```

### Similarity Heatmap

```python
spaxcr.pl.similarity_heatmap(
    cdr3s=cdr3sdf,
    label="Chain",
    similarity_method="gaussian",
    linkage_method="ward",
    insert_cost=1,
    replace_cost=1,
    triu=True,
)
```

### Radial Tree

```python
spaxcr.pl.similarity_radialtree(
    cdr3s=cdr3sdf["CDR3"],
    threshold=2,
)
```

Add metadata rings:

```python
spaxcr.pl.similarity_radialtree(
    cdr3s=cdr3sdf,
    threshold=2,
    cluster_in_meta=True,
    metas="Chain",
    insert_cost=2,
    replace_cost=1,
)
```

For large CDR3 sets, increase the number of worker processes:

```python
spaxcr.pl.similarity_heatmap(
    cdr3sdf,
    threads=4,
)
```

## 12. Sequence Logo

`spaxcr.pl.seqlogo` aligns protein sequences with ClustalW2 and plots a sequence logo with logomaker.

```python
clustalw_exe = "path/to/clustalw2"

spaxcr.pl.seqlogo(
    cdr3sdf["CDR3"],
    clustalw_path=clustalw_exe,
)
```

If `clustalw2` is available in your system `PATH`, `clustalw_path` can be omitted.

## 13. Query Clones and CDR3 Sequences

Find where a clone appears:

```python
spaxcr.tl.query(
    adata,
    var_names="clone_700",
)
```

Find where a CDR3 sequence appears:

```python
spaxcr.tl.query(
    adata,
    cdr3_names="CQTWGTGYVVF",
)
```

Find which clones are present in selected spots:

```python
spaxcr.tl.query(
    adata,
    obs_names=["AAACAAGTATCTCCCA", "AAACAGGGTCTATATT"],
)
```

Find all `IGH` clones in tumor spots:

```python
tumor_cells = adata[adata.obs["anno"] == "tumor"].obs_names

spaxcr.tl.query(
    adata,
    obs_names=tumor_cells,
    obs_chain="IGH",
    focus_cdr3=False,
)
```

Find all `IGHA` clones in tumor spots:

```python
spaxcr.tl.query(
    adata,
    obs_names=tumor_cells,
    obs_chain="IGHA",
    focus_cdr3=False,
)
```

## 14. Diversity Analysis

Compute diversity metrics by group:

```python
spaxcr.tl.diversity(
    adata,
    groupby="anno",
    chain="IGH",
    method="gini",
)
```

```python
spaxcr.tl.diversity(
    adata,
    groupby="anno",
    chain="IGH",
    method="shannon",
)
```

```python
spaxcr.tl.diversity(
    adata,
    groupby="anno",
    chain="IGH",
    method="simpson",
)
```

## 15. Density and Spatial Smoothing

Compute the density of a chain, clone set, or CDR3 set:

```python
spaxcr.tl.density(
    adata,
    chain="IGH",
    density_name="IGH_density",
    smooth_density=True,
    smooth_distance=500,
)
```

Smooth an existing numeric column in `adata.obs`:

```python
spaxcr.tl.smooth(
    adata,
    obs_name="sum_mu",
    alpha=1,
    kernel="Gaussian",
    smooth_distance=500,
    min_round_spot=1,
    fillna=False,
    add_obs="mu_density",
    norm_by="BCRclone_num",
)
```

Visualize the raw and smoothed values:

```python
sc.pl.spatial(
    adata,
    color=["sum_mu", "mu_density"],
    spot_size=200,
)
```

## 16. H/L or Alpha/Beta Spatial Pairing

If you have an alignment table with `x`, `y`, and `cloneId`, you can test spatial pairing between BCR heavy/light chains or TCR beta/alpha chains.

```python
pair_df = spaxcr.tl.XCR_pairing(
    adata,
    alignment=alignment_df,
    receptor_type="BCR",
    r=5,
    n_perm=1000,
    random_state=42,
    n_jobs=-1,
    min_observed=2,
)
```

For a general pairing test between two coordinate tables:

```python
pair_df = spaxcr.tl.spatial_pairing(
    df_a=df_a,
    df_b=df_b,
    r=5,
    n_perm=1000,
)
```

Both `df_a` and `df_b` must contain:

```text
name    x    y
```

## Next Steps

For a complete parameter reference, see [API Parameters](./api.md).
