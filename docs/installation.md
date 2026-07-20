# Installation

## Requirements

- Python >=3.10.0
- `bedtools` on PATH

## Quick install from GitHub

We recommend using a virtual environment to avoid dependency conflicts:

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade git+https://github.com/gopisiva1616/PyFuse.git
```


## Developer install from Git source

```bash
git clone https://github.com/gopisiva1616/PyFuse.git
cd PyFuse
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

Note: PyPI and Conda packages are not yet available. The commands below will be updated when the packages are published.

### Install from PyPI

```bash
pip install pyfuse
```

### Conda installation

```bash
conda install -c conda-forge pyfuse
```

## Verify installation

```bash
pyfuse -h
pyfuse annotator --help
pyfuse resources --help
```
