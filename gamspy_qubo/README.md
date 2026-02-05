# GAMSPy_QUBO: A Python Implementation of classic GAMS-QUBO reformulation tool

## Quickstart

### 1. Clone the repository

```bash
git clone git@git.gams.com:devel/qubo.git
cd qubo/gamspy_qubo
```

### 2. Create and activate a virtual environment

Note: This project uses `uv` as the environment management tool. It is recommended to use it with `uv`.

#### macOS / Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

#### Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
```

### 3. Install dependencies

#### Preferred method (using editable install):

```bash
pip install -e .
```

#### Fallback method (using `requirements.txt`):

```bash
pip install -r requirements.txt
```

## Project Structure

```
gamspy_qubo/
├── src/
│   └── gamspy_qubo/
│       ├── __init__.py
│       └── qubo.py
├── examples/
│   └── ...
├── tests/
│   └── ...
├── pyproject.toml
├── requirements.txt
├── README.md
└── .gitignore
```

## Using the Package

Once installed, you can use it like:

```python
from gamspy_qubo import Qubo
```

Or run example scripts:

```bash
python examples/qap.py
```

## Run Tests

There are tests available that can be run using the `pytest` package. You can install it via,

```bash
pip install -e .[test]
```
The tests can then be run via `pytest test_qubo.py`

## Build the Package (optional)

To create a distributable `.whl` or `.tar.gz`:

```bash
pip install build
python -m build
```

Artifacts will be available in the `dist/` directory.

## Classic GAMS-QUBO Reformulation tool

For information about the classic tool, please read [here](../README.md).