# SpaXCR

[![Documentation Status](https://readthedocs.org/projects/spaxcr/badge/?version=latest)](https://spaxcr.readthedocs.io/en/latest/)

SpaXCR is a toolkit for spatial single-cell TCR/BCR repertoire analysis.

## Documentation

The full documentation, including tutorials and API reference, is hosted on Read the Docs:

- [SpaXCR documentation](https://spaxcr.readthedocs.io/en/latest/)

## Installation

```bash
pip install spaxcr
```

## Quick start

```python
import spaxcr

adata = spaxcr.pp.load_xcr(
    clone_file="clone.txt",
    cloneId_barcode_file="clone_barcode.tsv",
)
```

The package exposes four main submodules:

- `spaxcr.pp`: preprocessing, loading, and filtering
- `spaxcr.tl`: analysis tools
- `spaxcr.pl`: visualization
- `spaxcr.extra`: additional helpers
