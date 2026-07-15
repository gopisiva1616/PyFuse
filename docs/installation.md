# Installation

## Requirements

- Python 3.10+
- `bedtools` on PATH

## Install from PyPI

```bash
pip install pyfuse
```

## Install from Git source

```bash
git clone https://github.com/gopisiva1616/PyFuse.git
cd PyFuse
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

## Conda installation

After conda package publication:

```bash
conda install -c conda-forge pyfuse
```

## Verify installation

```bash
pyfuse -h
pyfuse annotator --help
pyfuse resources --help
```
