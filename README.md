![Corromics banner](assets/corromics_full_logo.png)

# Corromics application

Corromics is an interactive Streamlit application for discovering and exporting
cross-omics associations. This repository contains the web application, Docker
configuration, native application build files, example data, and WSL/Linux
setup for optional Joint-RPCA/Gemelli analyses.

- Hosted application: <https://corromics.gnps2.org/>
- Application repository: <https://github.com/abzer005/Corromics>
- Python package: <https://github.com/abzer005/corromics_python_package>
- Windows installer: <https://www.functional-metabolomics.com/resources>

The application homepage contains the scientific workflow, input requirements,
method descriptions, and interpretation guidance. This README focuses on
installation, deployment, compatibility, and operations.

## Choose how to run Corromics

| Platform or use case | Recommended route | Alternative |
|---|---|---|
| No installation, small datasets | Hosted web app | Local Git clone |
| Windows, complete Linux-compatible workflow | WSL2 | Docker Desktop |
| Windows, standard workflow | Windows installer | Git clone, WSL2, or Docker Desktop |
| Linux | Git clone | Docker |
| macOS | Git clone | Docker Desktop |
| Notebooks, pipelines, servers, or HPC | Python package | Editable package clone |

## Hosted web application

For smaller datasets, open <https://corromics.gnps2.org/>. Hosted mode blocks
runs at 1,000,000 cross-omics associations. Use a local clone, Docker, or WSL
for larger analyses.

## Run from a Git clone on macOS, Linux, or Windows

Python 3.11 is recommended because it matches the application and Windows
build.

```bash
git clone https://github.com/abzer005/Corromics.git
cd Corromics
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run Home.py
```

Open <http://localhost:8501>. Stop the application with `Ctrl+C`. Closing the
browser does not stop the Streamlit process.

### Git not installed?

Check first:

```bash
git --version
```

If Git is missing:

- macOS: install Apple's command line tools with `xcode-select --install`, or install Git from <https://git-scm.com/download/mac>.
- Ubuntu/Debian Linux: run `sudo apt update && sudo apt install -y git`.
- Windows: install Git for Windows from <https://git-scm.com/download/win>, then reopen PowerShell and rerun `git --version`.

## Run with Docker

Local Docker users need these included files:

- `Dockerfile` builds the Corromics image.
- `run_server.sh` starts the application inside the container.
- `environment-gemelli-worker.yml` installs the Joint-RPCA/Gemelli environment
  during the image build.

The `docker-compose.yml` and `docker-compose-dev.yml` files are intended for the
maintainers of the hosted GNPS2 deployment. Local users do not need them.

See the [Docker Desktop user guide](Docker_User_Guide/README.md) for macOS and
Windows installation, build and run commands, Apple Silicon instructions,
memory guidance, upgrades, and troubleshooting.

## Run on Windows with WSL2

WSL2 is the recommended Windows route when Joint-RPCA/Gemelli is required.
Native Windows and WSL are separate installations.

The [`WSL_Joint_RPCA`](WSL_Joint_RPCA) folder contains everything needed for
this setup:

- [`README.md`](WSL_Joint_RPCA/README.md) provides the complete WSL instructions.
- [`setup_wsl_joint_rpca.sh`](WSL_Joint_RPCA/setup_wsl_joint_rpca.sh) installs
  the Corromics environment and Joint-RPCA/Gemelli dependencies.

Follow that folder's README and run the included shell script as shown below.

Open PowerShell as Administrator:

```powershell
wsl --install -d Ubuntu
```

Restart Windows if requested, open Ubuntu, and create the requested Linux user.
Then run:

```bash
sudo apt update
sudo apt upgrade -y
sudo apt install -y git
cd ~
git clone https://github.com/abzer005/Corromics.git
cd Corromics
bash WSL_Joint_RPCA/setup_wsl_joint_rpca.sh
conda activate corromics-main
streamlit run Home.py --server.port 5000 --server.address 127.0.0.1
```

Open <http://localhost:5000> in the Windows browser. Keep the repository under
the WSL home directory, such as `~/Corromics`, rather than `/mnt/c/...`; this is
usually faster and avoids cross-filesystem permission problems.

## Windows installer

Download the installer from the
[Functional Metabolomics Lab resources page](https://www.functional-metabolomics.com/resources),
run it, and launch Corromics from the installed shortcut.

The native Windows build supports the standard Corromics workflow. Joint-RPCA
is excluded because the full Gemelli dependency stack is not reliably available
as native Windows packages. Use WSL2 for Joint-RPCA.

## Python package for pipelines

The separate [`corromics_python_package`](https://github.com/abzer005/corromics_python_package)
repository provides the app-compatible analysis functions without Streamlit.
It is intended for Jupyter notebooks, repeatable batch processing, Snakemake or
Nextflow workflows, servers, and HPC environments.

Install it directly from GitHub:

```bash
pip install git+https://github.com/abzer005/corromics_python_package.git
```

Example:

```python
from corromics import pearson_omics, write_csv, write_graphml

results = pearson_omics(omics_1, omics_2, fdr=None)
write_csv(results, "pearson_results.csv")
write_graphml(results, "pearson_network.graphml")
```

The package repository documents its complete API, notebooks, optional
dependencies, testing, and output examples.

## Method availability

| Method or feature | Hosted web | Docker | WSL/Linux/macOS clone | Windows installer | Python package |
|---|---:|---:|---:|---:|---:|
| Hierarchical binning | Yes | Yes | Yes | Yes | Yes |
| Pearson | Yes | Yes | Yes | Yes | Yes |
| Spearman | Yes | Yes | Yes | Yes | Yes |
| Distance correlation | Yes | Yes | Yes | Yes | Yes |
| SparCC | Yes, subject to limits | Yes | Yes | Yes | Yes |
| Joint-RPCA | Yes | Yes | macOS/Linux/WSL with optional Gemelli worker | No | Linux/WSL with optional Gemelli environment |
| Target-decoy FDR | Pearson, Spearman, distance correlation | Same | Same | Same | Same helper functions |
| CSV and GraphML export | Yes | Yes | Yes | Yes | Yes |

SparCC and Joint-RPCA do not expose FDR filtering because the application does
not provide calibrated p-values for those methods.

## Hardware and memory recommendations

These are operational guidelines, not strict requirements:

| Workload | CPU | RAM | Notes |
|---|---:|---:|---|
| Hosted/small analysis | Browser only | Browser only | Server limits apply |
| Small local analysis | 2–4 cores | 8 GB | Suitable for exploratory runs |
| Medium local analysis | 4–8 cores | 16 GB | Recommended general configuration |
| Large local analysis | 8+ cores | 32 GB or more | Filter features before association analysis |
| SparCC or Joint-RPCA | 8+ cores | 32 GB or more | Memory use grows rapidly with features |

The supplied hosted Docker Compose configuration limits the service to 8 GB.
When using Docker Desktop, allocate at least 8 GB to Docker for normal runs and
more for large SparCC or Joint-RPCA analyses. The host operating system also
needs free memory.

## Large-dataset limits

The number of cross-omics pairs is approximately:

```text
number of omics-1 features × number of binned omics-2 features
```

- Hosted/restricted mode blocks runs at 1,000,000 associations.
- Local mode blocks runs at 10,000,000 associations to reduce crash risk.

For large analyses, use a local installation and reduce the feature count with
filtering or hierarchical binning if runtime or memory use becomes high.

## Troubleshooting

### Port already in use

Choose another Streamlit port:

```bash
streamlit run Home.py --server.port 8502
```

For Docker, map a different host port while keeping container port 5000:

```bash
docker run --rm -p 8502:5000 -e ENVIRONMENT_MODE=local corromics
```

Then open <http://localhost:8502>.

### Docker exits or is killed

- Increase Docker Desktop's memory allocation.
- Check available disk space with `docker system df`.
- Reduce the number of features before large SparCC analyses.
- Inspect a background container with `docker logs corromics`.
- Rebuild after dependency changes with `docker build --no-cache -t corromics .`.

### WSL problems

- Check WSL with `wsl --status` in PowerShell.
- Update it with `wsl --update`.
- Store the repository under `~/Corromics`, not `/mnt/c/...`.
- Rerun `bash WSL_Joint_RPCA/setup_wsl_joint_rpca.sh` if Gemelli is missing.
- Activate `corromics-main` before launching the application.
- If localhost forwarding fails, run `wsl --shutdown` in PowerShell, reopen
  Ubuntu, and launch Corromics again.

### Python installation errors

- Confirm the interpreter with `python --version`; Python 3.11 is recommended.
- Use `python -m pip` so dependencies install into the active interpreter.
- Use a fresh virtual environment instead of system Python.
- Upgrade tools with `python -m pip install --upgrade pip setuptools wheel`.
- Install `python3-venv` on Linux if virtual-environment creation fails.
- Native Windows Joint-RPCA dependency failures are expected; use WSL/Linux.

### `conda: command not found` in WSL/Linux

Check whether Miniconda or another Conda installation already exists:

```bash
ls -d ~/*conda*
```

If a Conda directory is listed, activate it using the path you found. For a
typical Miniconda install:

```bash
source ~/miniconda3/bin/activate
conda init bash
source ~/.bashrc
```

If Conda is not installed:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the installer prompts, initialize Conda when asked, then close and reopen
the shell or run `source ~/.bashrc`.

### Starting WSL after initial installation

Do not rerun the setup script for routine use. Open Ubuntu and run:

```bash
cd ~/Corromics
conda activate corromics-main
streamlit run Home.py --server.port 5000 --server.address 127.0.0.1
```

Then open <http://localhost:5000> in the Windows browser.

### Optional Joint-RPCA worker on macOS

Joint-RPCA can use a separate Python environment for Gemelli. After creating a
Gemelli-capable Conda environment, find its Python interpreter:

```bash
conda run -n gemelli-standalone python -c "import sys; print(sys.executable)"
```

Start Corromics with that interpreter path:

```bash
export GEMELLI_WORKER_PYTHON="/path/from/the/command/above"
export CORROMICS_GEMELLI_BACKEND=local
streamlit run Home.py
```

## Upgrading

Update a cloned application:

```bash
cd Corromics
git pull
source .venv/bin/activate
pip install -r requirements.txt
```

Update Docker:

```bash
cd Corromics
git pull
docker build --no-cache -t corromics .
```

Update the package:

```bash
pip install --upgrade --force-reinstall \
  git+https://github.com/abzer005/corromics_python_package.git
```

## Versioning and changelog

The application version is defined in `src/version.py`. This restructuring pass
prepares `v1.0.0` as the proposed first formal release, pending maintainer
review.

To reproduce an analysis later, record the Corromics version or Git commit
together with the selected method, settings, and random seed.

If the app and Python package are being compared, use versions released together
or record the exact commits used. This is the simplest way to obtain matching
results.

See [CHANGELOG.md](CHANGELOG.md) for release notes.

## Privacy statement for the hosted app

Files uploaded to the hosted app are processed on a remote Corromics/GNPS
server. If data confidentiality is a concern for your project, consider running
Corromics locally using a Git clone, WSL, Docker, or the Python package.

## Citation

Until a dedicated Corromics publication or archived software DOI is available,
cite the repository and exact version or commit used:

> Functional Metabolomics Lab. Corromics.
> <https://github.com/abzer005/Corromics>

Also cite the original publication for each analysis method used. Method
references are provided on the application homepage.

## License

Corromics is distributed under the BSD 3-Clause License. The license text is
available in [assets/corromics_license.rtf](assets/corromics_license.rtf).

## Contributors and contact

Corromics is developed by contributors at the
[Functional Metabolomics Lab](https://github.com/Functional-Metabolomics-Lab).

- Issues and feature requests: <https://github.com/abzer005/Corromics/issues>
- Lab website: <https://www.functional-metabolomics.com/>
- Contributions: open an issue or pull request in this repository.
