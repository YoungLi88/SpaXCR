from pathlib import Path

from setuptools import find_packages, setup


README = (Path(__file__).parent / "README.md").read_text(encoding="utf-8")


setup(
    name="spaxcr",
    version="0.1.0",
    description="SpaXCR: a toolkit for spatial single-cell TCR/BCR repertoire analysis",
    long_description=README,
    long_description_content_type="text/markdown",
    author="liuyihhha",
    author_email="liuyi212@mails.ucas.ac.cn",
    url="https://github.com/YoungLi88/SpaXCR",
    license="MIT",
    license_files=["LICENSE"],
    python_requires=">=3.8",
    install_requires=[
        "anndata",
        "scanpy",
        "pandas",
        "numpy",
        "scipy",
        "typing_extensions",
        "scikit-learn",
        "statsmodels",
        "tqdm",
        "matplotlib",
        "seaborn",
        "venn",
        "upsetplot",
        "logomaker",
        "biopython",
        "python-Levenshtein",
    ],
    packages=find_packages(include=["spaxcr", "spaxcr.*"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
