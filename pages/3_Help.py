import streamlit as st

from src.common import initialize_app
from src.version import __version__


st.set_page_config(
    page_title="CorrOmics Help",
    page_icon="❓",
    layout="wide",
)

initialize_app()


# ============================================================
# Header
# ============================================================

st.title("Help & Documentation")

st.markdown("""
Use this page for practical guidance on preparing input files, understanding
CorrOmics outputs, exporting results, running the app locally, and troubleshooting.

For complete installation and developer documentation, see the
[CorrOmics GitHub repository](https://github.com/abzer005/Corromics).
""")


# ============================================================
# Quick navigation
# ============================================================

st.subheader("Quick Guide")

st.markdown("""
| I want to... | Go to |
|---|---|
| Prepare my input files | [Input Data](#input-data) |
| Understand result columns | [Understanding Results](#understanding-results) |
| Understand target-decoy outputs | [Target-Decoy FDR](#target-decoy-fdr) |
| Export a network to Cytoscape | [GraphML & Cytoscape](#graphml-cytoscape) |
| Run CorrOmics locally | [Running CorrOmics](#running-corromics) |
| Set up Joint-RPCA | [Joint-RPCA Setup](#joint-rpca-setup) |
| Fix a problem | [Troubleshooting](#troubleshooting) |
""")


# ============================================================
# Input data
# ============================================================

st.markdown('<a id="input-data"></a>', unsafe_allow_html=True)
st.header("Input Data")

st.markdown("""
CorrOmics links metabolomics data with a complementary omics dataset using
a shared metadata table.

The application performs several formatting and matching steps automatically,
but sample identifiers must be consistent across the corresponding files.
""")

with st.expander("🧾 Metadata", expanded=False):
    st.markdown("""
The metadata table links samples between the metabolomics and complementary
omics datasets.

### Required columns

**`filename`**

Sample name corresponding to the metabolomics table.

Examples:

- `sample1.mzML`
- `sample2.mzML`

**`ATTRIBUTE_Corromics_filename`**

Corresponding sample identifier used in the complementary omics table.

Examples:

- sequencing filename
- microbiome sample ID
- proteomics sample ID

### Optional columns

Additional metadata can be included and used for sample filtering, for example:

- treatment
- control/sample type
- time point
- tissue
- condition
- batch

### Example

| filename | ATTRIBUTE_Corromics_filename | Sample_Type | Time_Point |
|---|---|---|---|
| sample1.mzML | sample1.fasta | Sample | 1h |
| sample2.mzML | sample2.fasta | Sample | 2h |
| sample3.mzML | sample3.fasta | Sample | 3h |
| blank.mzML | --- | Blank | N/A |

If samples fail to match, first check filename extensions, spaces,
capitalization, and prefixes/suffixes.
""")


with st.expander("🧬 Complementary Omics Feature Table", expanded=False):
    st.markdown("""
CorrOmics supports complementary quantitative omics tables such as:

- microbiome / ASV abundance
- metagenomics
- proteomics
- transcriptomics
- other feature-abundance datasets

### Expected structure

- features in **rows**
- samples in **columns**
- first column contains the feature identifier

Additional non-sample columns, such as taxonomy or annotations, may also be
included.

### Example

| feature_ID | sample1 | sample2 | sample3 | Genus |
|---|---:|---:|---:|---|
| ASV_1 | 120 | 340 | 220 | Bacteroides |
| ASV_2 | 80 | 20 | 150 | Blautia |

Sample columns must correspond to the identifiers provided in the metadata.
""")


with st.expander("⚛️ Metabolomics Feature Table", expanded=False):
    st.markdown("""
CorrOmics supports multiple metabolomics input formats.

### Manual Input

Expected format:

- features in **rows**
- samples in **columns**
- a column named `feature_ID`

Example:

| feature_ID | sample1.mzML | sample2.mzML | blank.mzML |
|---|---:|---:|---:|
| 1 | 1000 | 1100 | 100 |
| 2 | 2000 | 2200 | 200 |

Sample names should match the `filename` values in the metadata.

---

### FBMN-Stats Output

A supported FBMN-Stats feature table can be uploaded directly.

CorrOmics automatically restructures the supported output format for use in
the association workflow.

---

### MZmine Export

CorrOmics can process MZmine feature tables containing fields such as:

- `row ID`
- `row m/z`
- `row retention time`

These values are used to construct feature identifiers and identify the
sample-intensity columns.

---

### Example Dataset

CorrOmics also provides the **CCE1706 Metabolomics** example dataset for
exploring the workflow without preparing a metabolomics table manually.
""")


# ============================================================
# Understanding results
# ============================================================

st.markdown('<a id="understanding-results"></a>', unsafe_allow_html=True)
st.header("Understanding Results")

st.markdown("""
The exact result columns depend on the selected association method.
Not every method produces every statistical output.
""")


with st.expander("Result Table Columns", expanded=False):
    st.markdown("""
### `Feature`

Feature from one omics dataset.

### `Variable`

Feature from the paired omics dataset.

### `Estimate`

Association score produced by the selected method.

For Pearson and Spearman, this is a signed correlation coefficient.

For other methods, the interpretation depends on the method.

### `P-value`

Nominal statistical significance before multiple-testing correction.

### `BH-Corrected P-Value`

P-value adjusted using the Benjamini-Hochberg procedure.

### `R2`

Squared correlation coefficient where applicable.

### Important

Not all methods provide p-values, adjusted p-values, or R² values.

A missing value should not be interpreted as zero or as statistical
significance.
""")


with st.expander("Scientific Notation and Downloads", expanded=False):
    st.markdown("""
CorrOmics may display small numerical values using scientific notation to keep
large result tables readable.

The underlying numerical values are preserved.

Downloaded CSV files contain the numerical values used by CorrOmics rather
than shortened display strings.
""")


# ============================================================
# Target-decoy
# ============================================================

st.markdown('<a id="target-decoy-fdr"></a>', unsafe_allow_html=True)
st.header("Target-Decoy FDR")

with st.expander("What are target and decoy associations?", expanded=False):
    st.markdown("""
For supported association methods, CorrOmics can generate a decoy distribution
for empirical false discovery rate estimation.

Results are presented separately as:

- **Target associations**
- **Decoy associations**
- target association-score distribution
- decoy association-score distribution

The target table contains associations calculated from the original paired
data.

The decoy table provides a reference distribution generated from disrupted
feature relationships.

Comparing target and decoy score distributions allows CorrOmics to estimate
the background level of associations expected by chance.

Target-decoy FDR is not available for every association method.
""")


with st.expander("How should I use target-decoy results?", expanded=False):
    st.markdown("""
Target-decoy filtering is intended to help prioritize associations whose scores
are stronger than expected from the decoy background.

Users should still consider:

- effect size
- biological plausibility
- sample size
- preprocessing choices
- potential confounding
- independent validation

An association passing an FDR threshold does not demonstrate causality.
""")


# ============================================================
# Output files
# ============================================================

st.header("Output Files")


with st.expander("CSV Association Tables", expanded=False):
    st.markdown("""
CorrOmics association results can be downloaded as CSV files.

Depending on the analysis, downloads may include:

- complete target results
- complete decoy results
- filtered association tables
- method-specific result tables

These files can be used for downstream analysis in R, Python, Excel, or other
statistical software.
""")

st.markdown('<a id="graphml-cytoscape"></a>', unsafe_allow_html=True)
with st.expander("🌐 GraphML & Cytoscape", expanded=False):
    st.markdown("""
    ...
CorrOmics can export association networks in **GraphML** format for
visualization in Cytoscape.

### Importing into Cytoscape

1. Download the `.graphml` file from CorrOmics.
2. Open Cytoscape.
3. Select:

   `File → Import → Network from File…`

4. Select the downloaded GraphML file.

The network topology and available node/edge attributes should be imported
together.

### Node Table

CorrOmics may include node attributes such as:

- `Node_Info`
- `Original_index`

`Node_Info` can be used to distinguish the two omics blocks and assign
different node colors or shapes.

### Edge Table

Depending on the selected method, edge attributes may include:

- `Correlation_Score`
- `Sign_Score`
- `Absolute_Correlation_Score`
- `P_Value`
- `BH_Corrected_P_Value`
- `R2`
- `Method`

For signed methods, `Sign_Score` or `Correlation_Score` can be used to
distinguish positive and negative edges.

`Absolute_Correlation_Score` can be used to scale edge width.
""")


# ============================================================
# Interactive plots
# ============================================================

st.header("Interactive Plots")

with st.expander("Using Plotly Figures", expanded=False):
    st.markdown("""
Most CorrOmics figures are interactive.

You can:

- hover over points to inspect values
- drag to zoom into a region
- double-click to reset the view
- use the Plotly toolbar for additional controls
- export figures using the available download controls

If a calculation is still running, Streamlit may display a running indicator
in the application interface.
""")


# ============================================================
# Running CorrOmics
# ============================================================

st.markdown('<a id="running-corromics"></a>', unsafe_allow_html=True)
st.header("Running CorrOmics")

st.markdown("""
CorrOmics can be used through the hosted web application or installed locally.

The repository contains the complete and current installation instructions:

[github.com/abzer005/Corromics](https://github.com/abzer005/Corromics)
""")


with st.expander("🌐 Hosted CorrOmics", expanded=False):
    st.markdown("""
The hosted application is available at:

[https://corromics.gnps2.org/](https://corromics.gnps2.org/)

No installation is required.

The hosted application currently supports analyses with up to
**1,000,000 candidate associations**.

For larger analyses, use a local installation.
""")


with st.expander("💻 Local CorrOmics", expanded=False):
    st.markdown("""
Running CorrOmics locally provides a higher analysis limit and additional
deployment options.

Local CorrOmics currently supports up to:

**10,000,000 candidate associations**

Available installation routes include:

- **Windows installer**
- **Git-based installation**
- **Docker**
- **WSL2 for Windows**
- **Python package for programmatic workflows**

The recommended installation route depends on your operating system and
whether Joint-RPCA is required.

See the
[CorrOmics README](https://github.com/abzer005/Corromics)
for the current installation commands and requirements.
""")


with st.expander("🪟 Windows", expanded=False):
    st.markdown("""
### Windows installer

The standalone Windows installer is the simplest option for users who want to
run CorrOmics locally without managing a Python environment.

The installer can be obtained through the CorrOmics resources provided by the
Functional Metabolomics Lab.

### Joint-RPCA on Windows

Joint-RPCA is not currently supported through the native Windows executable.

Windows users who require Joint-RPCA should use:

- **WSL2**, or
- **Docker**

Detailed Windows and WSL2 instructions are available in the CorrOmics
repository.
""")


with st.expander("🍎 macOS / 🐧 Linux", expanded=False):
    st.markdown("""
CorrOmics can be run locally on macOS and Linux using the Git repository.

The repository contains instructions for:

- cloning CorrOmics
- creating the required Python environment
- installing dependencies
- starting the Streamlit application
- configuring the optional Gemelli worker for Joint-RPCA

See:

[CorrOmics GitHub repository](https://github.com/abzer005/Corromics)
""")


with st.expander("🐳 Docker", expanded=False):
    st.markdown("""
Docker provides a reproducible environment for running CorrOmics across
supported platforms.

This option can be useful when:

- you want an isolated installation
- local Python dependency management is inconvenient
- you need a reproducible deployment
- you require Joint-RPCA in an environment supported by the Docker setup

The repository contains a dedicated Docker guide with the current commands and
configuration.

See:

[CorrOmics GitHub repository](https://github.com/abzer005/Corromics)
""")


with st.expander("🐍 Python Package", expanded=False):
    st.markdown("""
CorrOmics also provides programmatic functionality for users who want to use
parts of the workflow from Python rather than through the Streamlit interface.

This option is intended primarily for:

- scripted analyses
- reproducible pipelines
- developers
- users integrating CorrOmics functionality into other workflows

Package installation and usage examples are maintained in the repository.
""")


# ============================================================
# Joint-RPCA
# ============================================================

st.markdown('<a id="joint-rpca-setup"></a>', unsafe_allow_html=True)
st.header("Joint-RPCA Setup")

with st.expander("Gemelli Worker Requirement", expanded=False):
    st.markdown("""
Joint-RPCA requires an optional **Gemelli worker environment**.

The main CorrOmics application and Gemelli may use separate Python
environments to avoid dependency conflicts.

### macOS / Linux

Follow the Gemelli worker installation instructions in the CorrOmics
repository.

The repository provides the environment definition required for the worker.

### Windows

Joint-RPCA is not currently available through the native Windows executable.

Use:

- **WSL2**, or
- **Docker**

for Joint-RPCA on Windows.

### More information

See the current installation instructions in:

[github.com/abzer005/Corromics](https://github.com/abzer005/Corromics)
""")


# ============================================================
# Limits
# ============================================================

st.header("Analysis Limits")

st.markdown("""
CorrOmics estimates the number of candidate pairwise associations from the
number of features in the two omics datasets.

For example:

`features in omics 1 × features in omics 2`

Current limits are:

| Environment | Candidate association limit |
|---|---:|
| Hosted CorrOmics | 1,000,000 |
| Local CorrOmics | 10,000,000 |

If the analysis exceeds the available limit, reduce the feature space using
appropriate filtering before running the association analysis.
""")


# ============================================================
# Troubleshooting
# ============================================================

st.markdown('<a id="troubleshooting"></a>', unsafe_allow_html=True)
st.header("Troubleshooting")


with st.expander("Samples are not matching", expanded=False):
    st.markdown("""
Check that sample identifiers are consistent across:

- metadata
- metabolomics table
- complementary omics table

Common causes include:

- different filename extensions
- leading or trailing spaces
- differences in capitalization
- extra prefixes or suffixes
- different sample naming conventions
- missing samples in one of the datasets

The metadata must provide the correct correspondence between the two omics
tables.
""")


with st.expander("My analysis exceeds the association limit", expanded=False):
    st.markdown("""
The number of candidate pairwise associations increases rapidly with feature
count.

For example:

- 1,000 metabolites × 500 microbial features = 500,000 associations
- 5,000 metabolites × 2,000 microbial features = 10,000,000 associations

Consider biologically appropriate filtering such as:

- prevalence filtering
- abundance/intensity filtering
- metadata-based sample selection
- feature annotation-based filtering

Avoid filtering solely to obtain a desired statistical result.
""")


with st.expander("Joint-RPCA is unavailable", expanded=False):
    st.markdown("""
Joint-RPCA requires the optional Gemelli worker environment.

If CorrOmics cannot detect Gemelli:

1. confirm that the Gemelli worker environment has been installed
2. restart CorrOmics after installation
3. check the operating-system-specific instructions in the repository

For Windows, use WSL2 or Docker when Joint-RPCA is required.
""")


with st.expander("The app appears to be stuck or still calculating", expanded=False):
    st.markdown("""
Large association calculations may take time depending on:

- number of candidate associations
- selected association method
- computer memory
- available CPU resources

Check whether the Streamlit interface indicates that the application is still
running.

If the calculation repeatedly fails, reduce the analysis size or consult the
repository troubleshooting documentation.
""")


with st.expander("Something else is not working", expanded=False):
    st.markdown("""
First check the latest documentation in the CorrOmics repository:

[github.com/abzer005/Corromics](https://github.com/abzer005/Corromics)

If the problem persists, report it through the GitHub Issues page.

When reporting a problem, include where possible:

- CorrOmics version
- operating system
- installation method
- selected association method
- error message
- steps required to reproduce the problem

Do not upload confidential or sensitive datasets to a public GitHub issue.
""")


# ============================================================
# Documentation
# ============================================================

st.header("More Documentation")

st.markdown("""
The **CorrOmics GitHub repository** is the main source for detailed technical
documentation.

Use the repository when you need:

- complete installation commands
- Windows installation information
- macOS or Linux setup
- Docker setup
- WSL2 setup
- Joint-RPCA / Gemelli configuration
- updating an existing installation
- Python package information
- release and version information
- developer documentation

**Repository:**  
[github.com/abzer005/Corromics](https://github.com/abzer005/Corromics)

This Help page is intended as a practical guide for using CorrOmics rather than
a duplicate of the full repository documentation.
""")


# ============================================================
# Feedback
# ============================================================

st.header("Feedback & Issues")

st.markdown("""
Feedback, bug reports, and feature requests can be submitted through the
CorrOmics GitHub repository.

[Open CorrOmics on GitHub](https://github.com/abzer005/Corromics)

For technical problems, please provide enough information to reproduce the
issue where possible.
""")


# ============================================================
# Footer
# ============================================================

st.markdown("---")
st.caption(f"CorrOmics v{__version__}")