---
title: Installation
description: Install Dense Arrays in a Python environment, add optional features, and configure FIMO when a preparation recipe needs it.
---

# Install Dense Arrays

Use Python 3.12 or later. The base package includes NumPy, OR-Tools with the
CBC packing backend, and the `dense-arrays` command. Start with pip in a
virtual environment; optional project-managed environments are described below.

## Install the published package

From a directory where you keep analysis projects, run these commands in a
macOS or other POSIX shell:

```bash
# Create a directory for your analysis and enter it.
mkdir array-project
cd array-project

# Create an isolated environment with Python 3.12.
python3.12 -m venv .venv

# Select this environment's Python and commands for the current shell.
source .venv/bin/activate

# Install Dense Arrays and its required Python dependencies from PyPI.
python -m pip install dense-arrays

# Show the installed command's options.
dense-arrays --help
```

If your compatible Python has another executable name, replace `python3.12`
with that name. Use `source .venv/bin/activate` again in each new terminal.
The [Python Packaging guide](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/)
explains virtual environments. On Windows PowerShell, activate the environment
with `.venv\Scripts\Activate.ps1` instead of the `source` command.

Continue with [your first array](quickstart.md). That example uses the published
optimizer and needs no FIMO installation.

The published [0.2.1 package](https://pypi.org/project/dense-arrays/0.2.1/)
provides direct packing and playback. The parts-to-library workflow and `tables`
extra described in this documentation require the source checkout containing
those features; they are not in the published 0.2.1 distribution.

## Optional features

Choose the dependencies required by your task:

| Task | Installation | Availability |
| --- | --- | --- |
| Render placement images and GIFs | `python -m pip install "dense-arrays[playback]"` | Published package |
| Write MP4 playback | Install `playback` and an FFmpeg executable on `PATH` | Published package |
| Read Parquet and Excel part tables | `python -m pip install ".[tables]"` from the source checkout | Source workflow |
| Use CSV or TSV part tables | No table extra | Source workflow |
| Score sampled candidates with FIMO | Install the MEME Suite executable separately | Source workflow |

The `docs` and `dev` extras support documentation and contributor checks; see
[development](development.md). FIMO is an external executable, not an extra.

## Use the library workflow

### Install from source

The library workflow is available on `feat/library-workflow`. Clone that branch,
then create an environment for the checkout. These commands require Git:

```bash
# Download the source that provides the library workflow.
git clone --branch feat/library-workflow https://github.com/e-south/dense-arrays.git
cd dense-arrays

# Create and select an environment for this checkout.
python3.12 -m venv .venv
source .venv/bin/activate

# Install the checkout with optional table readers and playback rendering.
python -m pip install ".[tables,playback]"

# Check that this installation exposes workflow planning.
dense-arrays plan --help
```

Use `python -m pip install .` if you do not need those extras. This installs a
snapshot of the local source; rerun the installation after updating the checkout.
If you already have the workflow checkout, start with the environment commands
from its root directory.

## Configure FIMO for motif scoring

Preparation recipes that request `FimoScoring` need the MEME Suite `fimo`
command. Dense Arrays resolves it from `PATH`, or accepts an explicit executable
path. Importing the Python package does not invoke FIMO.

Install the command-line tools using the
[MEME Suite installation guide](https://meme-suite.org/meme/doc/install.html).
The [Pixi route below](#keep-python-and-fimo-together-with-pixi) keeps FIMO and
Python in one project environment. If FIMO is already installed, check it from
the shell where you will run Dense Arrays:

```bash
fimo --version  # Confirm that the configured executable can run.
```

Follow the installer's `PATH` instructions if the shell cannot find `fimo`.
Alternatively, set the executable location in Python; replace the example path
with the installed binary:

```python
from dense_arrays import parts

# Bind subsequent preparation planning to this particular executable.
scoring = parts.FimoScoring(executable="/absolute/path/to/fimo")
```

Planning checks the executable and records its version and fingerprint.
Scoring uses FIMO's text output, configured threshold, background and strand
policy. See [motif scoring](reference/motif-scoring.md) for those settings.

`pymemesuite` supplies a [Python interface to MEME internals](https://pypi.org/project/pymemesuite/).
It does not satisfy this executable-based scoring interface. Installing it is
not a substitute for installing `fimo`.

## Manage an analysis project with uv

Use [uv projects](https://docs.astral.sh/uv/guides/projects/) when you want a
dependency declaration and lockfile alongside your analysis scripts. Start from
the directory that will contain the new project:

```bash
# Create an independent analysis project targeting Python 3.12.
uv init --no-package --no-workspace --python 3.12 array-project
cd array-project

# Add the published package and optional rendering dependencies.
uv add "dense-arrays[playback]"

# Run the installed command in the project's environment.
uv run dense-arrays --help
```

`pyproject.toml` declares the dependencies; `uv.lock` records their resolved
versions. uv creates `.venv` as needed and selects it for `uv run`, so shell
activation is unnecessary. Keep the manifest and lockfile with your analysis.
The `--no-workspace` option keeps this project separate from any parent project.
Run an existing Python script with `uv run python path/to/analysis.py`.

To use a local workflow checkout instead of the published package, replace the
`uv add` command with this [path dependency](https://docs.astral.sh/uv/concepts/projects/dependencies/#path):

```bash
# Replace ../dense-arrays with the path to the checkout containing the workflow.
uv add "dense-arrays[tables,playback] @ ../dense-arrays"

# Confirm that the selected source provides workflow planning.
uv run dense-arrays plan --help
```

## Keep Python and FIMO together with Pixi

Use a [Pixi project](https://pixi.prefix.dev/latest/python/tutorial/) if you want
one environment for Python packages and the MEME Suite executable. Bioconda's
[`meme` package](https://bioconda.github.io/recipes/meme/README.html) includes FIMO.
Choose a new project directory. Replace the absolute file URL below with the
location of a checkout containing the workflow:

```bash
# Create a project using the channels that supply Python and MEME Suite.
pixi init --channel conda-forge --channel bioconda array-project
cd array-project

# Resolve Python and the MEME Suite command-line programs together.
pixi add "python=3.12" meme

# Add the local Python package through Pixi's PyPI dependency support.
pixi add --pypi "dense-arrays @ file:///absolute/path/to/dense-arrays"

# Confirm both commands are available in the same project environment.
pixi run fimo --version
pixi run dense-arrays plan --help
```

Pixi records dependencies in `pixi.toml` and resolved packages in `pixi.lock`.
Run analysis commands through `pixi run` so they see the environment's FIMO.
The commands follow the official [initialization](https://pixi.prefix.dev/latest/reference/cli/pixi/init/)
and [dependency](https://pixi.prefix.dev/latest/reference/cli/pixi/add/) interfaces.
Keep the manifest and lockfile with your analysis.
