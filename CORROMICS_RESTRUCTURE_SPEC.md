# CORROMICS_RESTRUCTURE_SPEC.md

## 1. Goal

Restructure and clean up CorrOmics without changing its core scientific idea.

This pass should:
- improve warnings and method guidance;
- correct SparCC behavior/documentation for the current cross-omics workflow;
- make Joint-RPCA setup/errors operating-system aware;
- correct hosted-vs-local deployment detection and pair limits;
- improve result-table formatting and selected UI layout;
- make installation/documentation easier to navigate;
- audit repository structure and redundancy;
- preserve compatibility with the GNPS Streamlit/Docker template;
- prepare CorrOmics for formal versioned releases.

Do not redesign the app into a different product.

---

## 2. Required pre-edit audit

Before changing code:

1. Read the full `README.md`.
2. Inspect every repository folder and important root-level file.
3. Categorize files into:
   - app runtime;
   - analysis/scientific code;
   - third-party/vendor;
   - hosted GNPS/Docker deployment;
   - local Docker;
   - Windows/WSL;
   - Windows packaging/build;
   - development/configuration;
   - docs;
   - example data/assets;
   - generated/temporary;
   - candidate redundant/obsolete.
4. Compare the repo with:
   `https://github.com/Wang-Bioinformatics-Lab/Template_Streamlit_App`
5. Preserve the template-critical deployment files and conventions required by GNPS.
6. Report proposed removals before deleting anything.

The current CorrOmics root already contains template-like files including Dockerfile, Makefile, docker-compose files, `run_server.sh`, `pages/`, and requirements. Do not simplify these blindly if GNPS relies on them.

---

## 3. SparCC

### 3.1 Scientific behavior

The current cross-omics use case combines metabolomics and a complementary omics table, typically microbiome data.

For the intended metabolite–microbiome workflow:
- metabolomics may already be zero-imputed and TIC-normalized;
- microbiome may be converted to relative abundance;
- a SparCC-compatible relative-abundance option should handle zeros before relative-abundance normalization;
- the combined matrix is passed to the vendored SparCC implementation.

Use:

```python
corr_values, _ = sparcc3_module.main(
    counts_samples_by_features,
    method="SparCC",
    iter=1,
    norm="normalize",
    xiter=10,
    oprint=False,
)
```

`xiter=10` may remain implicit if that is the vendored default.

### 3.2 Remove outer iteration UI

Remove the user-facing control for "SparCC iterations."

The outer iteration count must be fixed at `1` because `norm="normalize"` is deterministic. Do not show text claiming that additional outer iterations improve stability.

Do not remove SparCC's internal sparse exclusion iterations.

### 3.3 Warning/guidance text

Remove obsolete warning text such as:

> SparCC is intended for compositional datasets ... It is not recommended for general metabolomics or proteomics intensity tables ...

and iteration-focused messages such as:

> More iterations can make the SparCC estimate more stable...

Replace with a concise, accurate method note along the following lines:

> **SparCC is a compositional log-ratio association method originally developed for within-microbiome abundance data. In CorrOmics, its use for metabolite–microbiome associations is exploratory because the two omics blocks originate from different measurement processes. Use non-negative, appropriately normalized inputs and interpret cross-omics SparCC associations cautiously.**

Do not claim that SparCC performs CLR internally.

### 3.4 Memory warning

Keep memory protection if useful, but remove iteration-related advice.

The memory warning should focus on:
- number of combined features;
- estimated memory risk;
- available memory;
- filtering features before running large analyses.

Do not say "reduce iterations" because the outer iteration count is fixed at 1.

Review whether the current memory estimate constants still make sense when `iter=1` is fixed.

### 3.5 Preprocessing option

Keep the existing:
- `Relative abundance`: zero-preserving relative abundance.

Add or retain a distinct SparCC-compatible option such as:
- `Relative abundance (+1 pseudocount)`

Behavior for rows=features, columns=samples:

```python
X = df_transformed[cols_to_transform].astype(float) + 1
col_sums = X.sum(axis=0).replace(0, np.nan)
rel = X.div(col_sums, axis=1)
```

Do not silently change existing zero-preserving relative abundance.

Update the transformation caption accordingly.

---

## 4. Joint-RPCA / Gemelli

### 4.1 Worker resolution

The app must support a separate Gemelli Python environment.

Resolution order:
1. If backend is disabled, return no worker.
2. If `GEMELLI_WORKER_PYTHON` is set:
   - resolve/expand path;
   - ensure it exists;
   - verify the interpreter can import `gemelli`;
   - use it to launch `src/joint_rpca_worker.py`.
3. Otherwise, if the currently running Python can import Gemelli, use `sys.executable`.
4. Otherwise return a clear dependency failure.

The existing Docker environment variable `GEMELLI_WORKER_PYTHON` must actually be honored.

### 4.2 OS-aware errors

macOS:
> Joint-RPCA requires the optional Gemelli worker environment. Create/configure the Gemelli environment and point CorrOmics to its Python interpreter. See the installation instructions.

Do not mention WSL as a requirement on macOS.

Linux/WSL:
> Joint-RPCA requires the optional Gemelli worker environment. Install/configure Gemelli for this CorrOmics installation.

Native Windows:
> Joint-RPCA is not supported by the native Windows installation. Use WSL2 or Docker Desktop.

Docker/hosted:
If the worker is missing, report a deployment/configuration problem rather than telling the user to install WSL.

Avoid duplicating OS-detection logic across several files.

---

## 5. Hosted vs local mode and association limits

Required limits:
- hosted/restricted GNPS mode: **1,000,000** cross-omics pairs;
- local mode: **10,000,000** cross-omics pairs.

Current problem:
The GNPS app runs in Docker but can appear as "local mode."

Fix this by using an explicit deployment-mode signal rather than simply assuming that Docker means local.

Recommended design:
- one environment variable, e.g. `ENVIRONMENT_MODE` or an existing equivalent;
- allowed values should clearly distinguish `hosted`/`restricted` and `local`;
- hosted production compose/config must explicitly set hosted/restricted;
- local Docker instructions must explicitly set local where needed;
- local Git clone should default to local;
- no ambiguous sidebar label.

Sidebar should clearly state the active limit, for example:
- `Hosted mode · limit: 1,000,000 associations`
- `Local mode · limit: 10,000,000 associations`

Review all places calculating or displaying `max_corr` so the underlying restriction and displayed text cannot disagree.

---

## 6. Association result tables

Apply consistent numeric display formatting across association-result tables.

### Display
Use scientific notation for decimal numerical values, especially:
- Estimate where appropriate;
- P-value;
- BH-corrected P-value;
- R2 where applicable;
- other continuous score columns.

Avoid converting the underlying numeric dataframe columns into formatted strings if doing so would break sorting/filtering/downloads.

Use dataframe display configuration/column formatting where possible.

### Full p-value inspection
Requirement:
- normal table display: scientific notation;
- when the user inspects/edits/double-clicks a p-value cell, the underlying full numeric precision should remain available.

Implement this using Streamlit-supported table configuration if feasible.
If Streamlit's current dataframe widget cannot guarantee a special double-click representation, preserve the numeric value and provide the closest reliable inspection behavior rather than faking a feature. Report any UI limitation.

CSV exports must retain numeric values, not formatted scientific-notation strings unless the existing export intentionally formats them.

---

## 7. UI restructuring

### 7.1 Target and decoy dataframes
Place:
- Target dataframe
- Decoy dataframe

in tabs rather than vertically stacked sections.

### 7.2 Target and decoy score tables
Place:
- association/correlation scores of target dataframe
- association/correlation scores of decoy dataframe

in tabs.

Keep downloads associated with the appropriate tab/table.

### 7.3 Example dataset label
Rename:

`Use Example Dataset`

to:

`Example Dataset - CCE1706 Metabolomics`

Check all duplicated labels/help text so the same naming is used consistently.

---

## 8. Terminology

Use **association** as the general umbrella term.

Recommended:
- "Association Analysis"
- "Association methods"
- "Association results"
- "Number of associations"

Use **correlation** when mathematically specific:
- Pearson correlation
- Spearman correlation
- SparCC correlation coefficient
- correlation matrix
- target/decoy correlation score when that object is specifically a correlation coefficient.

Joint-RPCA should not be generically described as a conventional correlation coefficient if that is not what the output represents.

Do not perform a blind global string replacement from "correlation" to "association."

---

## 9. Homepage

### 9.1 Method count
Change:

> CorrOmics currently supports four association methods

to:

> CorrOmics currently supports five association methods

Verify that all five methods are actually listed:
- Pearson
- Spearman
- Distance correlation
- SparCC
- Joint-RPCA

### 9.2 Cytoscape / GraphML guidance

Retain/expand:
`Tips for Using the GraphML file in Cytoscape`

Explain the direct GraphML route:
1. Download GraphML from CorrOmics.
2. In Cytoscape: `File → Import → Network from File…`
3. Select the `.graphml` file.
4. Cytoscape should import network topology and GraphML attributes.

Also explain:
- where the **Node Table** is visible in Cytoscape;
- where the **Edge Table** is visible;
- which CorrOmics attributes map to node attributes and edge attributes;
- how users can use edge attributes such as association estimate, p-value/FDR where available, method, source/target type, etc. in Cytoscape Styles.

Before writing attribute-specific documentation, inspect the actual CorrOmics GraphML writer/export structure and document only attributes that are really exported.

Also mention CSV import as an alternative:
`File → Import → Network from File…`, then map source and target columns and retain additional columns as edge attributes.

---

## 10. README and installation documentation

README should remain the primary place for installation/deployment/operations.

Homepage should provide concise user-facing guidance and links/references, not duplicate the entire README.

### 10.1 Local Linux/macOS clone
Use:

```bash
git clone https://github.com/abzer005/Corromics.git
cd Corromics
```

Keep a clear supported Python recommendation consistent with actual dependency constraints.

### 10.2 WSL setup

Document initial setup:

```bash
cd ~
git clone https://github.com/abzer005/Corromics.git
cd Corromics
bash WSL_Joint_RPCA/setup_wsl_joint_rpca.sh
conda activate corromics-main
streamlit run Home.py --server.port 5000 --server.address 127.0.0.1
```

### 10.3 `conda: command not found`

Add troubleshooting.

First:

```bash
ls -d ~/*conda*
```

If Miniconda is already installed, explain the appropriate activation/init sequence, e.g.:

```bash
source ~/miniconda3/bin/activate
conda init bash
source ~/.bashrc
```

If it is not installed:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Then follow installer prompts, initialize Conda, and reload the shell.

Do not assume the installation path is always `~/miniconda3`; mention that users should use the path found by the check if different.

### 10.4 Starting WSL after initial installation

Add a short repeat-use section such as:

```bash
cd ~/Corromics
conda activate corromics-main
streamlit run Home.py --server.port 5000 --server.address 127.0.0.1
```

Then open `http://localhost:5000`.

Do not make users rerun the setup script every time.

### 10.5 Joint-RPCA on macOS
Document the optional separate Gemelli worker environment and `GEMELLI_WORKER_PYTHON`.

Do not hard-code a user-specific `/opt/homebrew/...` path.

Provide a dynamic resolution command such as:

```bash
conda run -n gemelli-standalone python -c "import sys; print(sys.executable)"
```

### 10.6 Docker and GNPS
Keep GNPS template files intact unless comparison with the upstream template proves a safe improvement.

Clearly separate:
- hosted GNPS production Docker;
- local Docker Desktop;
- Docker development/maintainer workflow.

---

## 11. GNPS template compatibility

Compare CorrOmics against:
`Wang-Bioinformatics-Lab/Template_Streamlit_App`

Template-relevant files currently include:
- `Dockerfile`
- `Makefile`
- `docker-compose.yml`
- `docker-compose-dev.yml`
- `requirements.txt`
- `run_server.sh`
- app homepage entry point
- `pages/`

Do not rename/remove these merely for stylistic cleanup.

Where CorrOmics has additional files required for:
- Gemelli/Joint-RPCA;
- Windows packaging;
- local installation;
- app-specific dependencies;

retain them if functional.

Document intentional differences from the template.

---

## 12. Repository redundancy audit

Inspect all current root folders/files, including at minimum:

- `.devcontainer/`
- `.github/workflows/`
- `.streamlit/`
- `Docker_User_Guide/`
- `WSL_Joint_RPCA/`
- `assets/`
- `example-data/`
- `hooks/`
- `pages/`
- `src/`
- `.gitignore`
- `Dockerfile`
- `Home.py`
- `Makefile`
- `README.md`
- `README_JOINT_RPCA.txt`
- `docker-compose-dev.yml`
- `docker-compose.yml`
- `environment-gemelli-worker.yml`
- `environment.lock.yml`
- `install_corromics.sh`
- `pyproject.toml`
- `requirements.lock.txt`
- `requirements.txt`
- `run.py`
- `run_app_temp.spec`
- `run_server.sh`

For each candidate redundancy, determine:
- whether it is referenced by code;
- whether it is referenced by docs;
- whether it is used in CI;
- whether it is required by GNPS;
- whether it is required for Windows packaging;
- whether it is required for WSL/Joint-RPCA;
- whether it is generated and should instead be ignored.

Do not delete anything until the user approves the proposed-removal list.

---

## 13. Licensing

Current README states BSD 3-Clause as the CorrOmics license.

Do not finalize license restructuring in this pass without user approval.

Perform an audit:
1. Identify CorrOmics-owned code.
2. Identify vendored/third-party code.
3. Find license/copyright notices for each vendored component.
4. Determine which notices must remain.
5. Identify any conflicting or unclear top-level license files.
6. Recommend a clean final layout:
   - top-level `LICENSE` for CorrOmics;
   - third-party notices/licenses retained near vendor code or in an appropriate notices location.

Never replace a third-party license simply because the app uses BSD 3-Clause.

---

## 14. Version control and releases

### Current state to account for
The public GitHub repository currently shows no formal Releases/tags.
The current `pyproject.toml` contains build-system configuration but no application version.

### Introduce a release-ready versioning system

Use Semantic Versioning unless there is a strong repo-specific reason not to.

Recommended first formal release after this cleanup:
- propose `v1.0.0` if this is intended to represent the first stable public manuscript/app release;
- otherwise propose `v0.x.0` if the user considers the app pre-stable.

Do not decide silently. Report the recommendation before tagging.

Use one source of truth for the app version, for example:
- `src/version.py`, or
- another small dedicated version module compatible with the current non-package Streamlit layout.

Expose the version somewhere low-noise, such as:
- homepage footer/about area, and/or
- sidebar footer.

Add a `CHANGELOG.md` if absent, using a simple structure:
- Added
- Changed
- Fixed
- Documentation
- Known limitations

For this restructuring, release notes should include at least:
- SparCC cross-omics handling/documentation changes;
- removal of user-facing SparCC outer iterations;
- OS-aware Joint-RPCA worker guidance;
- hosted/local mode detection fix;
- association table formatting;
- target/decoy tab layout;
- Cytoscape documentation expansion;
- installation/README cleanup;
- repository/deployment cleanup if performed.

Do not automatically:
- create Git tag;
- create GitHub Release;
- push;
- publish installer artifacts.

At the end provide the exact commands the user can run to create the release after review.

---

## 15. Local-only Codex instruction/spec files

These planning files do not need to be committed.

Preferred local files:
- `AGENTS.md`
- `CORROMICS_RESTRUCTURE_SPEC.md`

If the user wants them invisible to Git without changing the shared `.gitignore`, add them to the repository's local:

`.git/info/exclude`

Example entries:

```text
AGENTS.md
CORROMICS_RESTRUCTURE_SPEC.md
```

This keeps the files local to that clone and does not alter the repository's tracked `.gitignore`.

If either file is already tracked, `.git/info/exclude` will not untrack it; report that case before proceeding.

---

## 16. Completion criteria

Before declaring completion:

### Scientific
- SparCC uses `norm="normalize"`, outer `iter=1`.
- SparCC iteration UI removed.
- SparCC warning updated accurately.
- Internal SparCC exclusion algorithm retained.
- Joint-RPCA science unchanged.
- OS-aware Joint-RPCA setup/errors work.

### UI
- target/decoy tables use tabs.
- target/decoy score tables use tabs.
- example dataset renamed.
- homepage says five methods.
- generic workflow terminology uses "association" appropriately.
- result tables use scientific-notation display without destroying numeric data.

### Deployment
- hosted GNPS = 1M limit.
- local = 10M limit.
- hosted Docker no longer incorrectly reports local mode.
- template-critical GNPS files preserved.

### Docs
- README installation routes are clear.
- WSL first-install and repeat-start instructions are distinct.
- Conda-not-found guidance added.
- macOS Joint-RPCA guidance corrected.
- Cytoscape GraphML section expanded from actual exported attributes.

### Repository
- complete redundancy report produced.
- no unapproved destructive cleanup.
- licensing audit produced.
- versioning/changelog prepared.

### Final report
Provide:
- changed files;
- proposed removals;
- retained files and reasons;
- tests;
- unresolved issues;
- proposed release version;
- release commands;
- `git diff --stat`;
- key diffs requiring review.


---

## 17. Simplify installation choices for non-technical scientific users

The installation/help structure should be simplified around **how a scientist wants to access CorrOmics**, rather than exposing many technical deployment details up front.

### 17.1 Present four clear access paths

Use a compact overview such as:

1. **Web app — easiest, no installation**
   - Works from any modern operating system/browser.
   - Best for users who want to run CorrOmics without local setup.
   - Clearly state the hosted service limitation: maximum **1,000,000 associations**.
   - Mention that large analyses may need local installation.

2. **Git clone — recommended local installation for macOS, Linux, and Windows**
   - Git itself is cross-platform; cloning the repository is valid on macOS, Linux, and Windows provided Git is installed.
   - Present this as the main local-development / local-scientist route.
   - Local mode allows up to **10,000,000 associations**, subject to the user's machine resources.
   - Explain OS-specific launch/setup only after the common clone step.

   Common start:
   ```bash
   git clone https://github.com/abzer005/Corromics.git
   cd Corromics
   ```

   Then branch into short OS-specific instructions:
   - **macOS**: local Python/Conda environment; optional separate Gemelli worker for Joint-RPCA.
   - **Linux**: local Python/Conda environment; optional Gemelli worker for Joint-RPCA.
   - **Windows**: repository can be cloned natively, but clearly distinguish which methods are supported natively and which require WSL2/Docker.

   Do not imply that Git clone itself is a macOS/Linux-only feature.

3. **WSL2 — Windows route when Joint-RPCA is needed**
   - Present WSL2 as a **Windows-specific compatibility route**, not as a general installation method all users need.
   - Explain that native Windows users who need Joint-RPCA should use WSL2 or Docker.
   - Reuse the existing `WSL_Joint_RPCA/` setup script and documented workflow.
   - Avoid making macOS/Linux users read WSL instructions.

4. **Docker — portable containerized route**
   - Explain that Docker can be used on macOS, Linux, and Windows.
   - Do not state that Docker is only or primarily a Windows solution.
   - It may be particularly convenient on Windows when users want an isolated environment and/or to avoid native dependency issues, but Docker Desktop's Windows backend details should not be presented as something scientists need to understand.
   - Keep the explanation user-oriented: "use Docker if you prefer an isolated, reproducible installation."
   - Distinguish local Docker from the hosted GNPS Docker deployment.

### 17.2 Use a simple decision guide

Add a compact decision guide near the installation section, for example:

| User need | Recommended route |
|---|---|
| I just want to use CorrOmics | Web app |
| I want to run larger analyses locally | Git clone |
| I use macOS/Linux and need all locally supported methods | Git clone + optional Gemelli worker |
| I use Windows and need Joint-RPCA | WSL2 or Docker |
| I want an isolated/reproducible environment | Docker |

Do not overwhelm users with implementation details in this overview.

### 17.3 Progressive disclosure

The README/Home page should use progressive disclosure:

- First show **which option should I choose?**
- Then show **minimal commands**.
- Put troubleshooting and advanced technical details below or in expandable/clearly separated sections.

Target audience:
- scientists comfortable with analysis tools;
- not necessarily comfortable with Conda internals, Docker architecture, Linux shells, environment variables, or deployment terminology.

Avoid unexplained terms such as:
- "backend"
- "container runtime"
- "WSL distribution"
- "environment resolver"

when a simpler phrase works.

### 17.4 Method-availability map

Provide a concise method-availability table tied to access route.

At minimum verify and document:
- Pearson
- Spearman
- Distance correlation
- SparCC
- Joint-RPCA

The table must be based on actual current implementation and tested environment behavior, not assumptions.

Suggested structure:

| Access route | Pearson | Spearman | Distance correlation | SparCC | Joint-RPCA | Association limit |
|---|---:|---:|---:|---:|---:|---:|
| Web app | ✓ | ✓ | ✓ | ✓ | verify | 1M |
| macOS Git clone | ✓ | ✓ | ✓ | ✓ | ✓ with Gemelli worker | 10M |
| Linux Git clone | ✓ | ✓ | ✓ | ✓ | ✓ with Gemelli worker | 10M |
| Windows native Git clone | verify | verify | verify | verify | ✗ | 10M if local |
| Windows WSL2 | ✓ | ✓ | ✓ | ✓ | ✓ | 10M |
| Local Docker | verify/test | verify/test | verify/test | verify/test | verify/test | 10M |

Do not publish checkmarks until verified against the repository and tested/known runtime behavior.

### 17.5 Home page wording

The Home page should give the shortest usable orientation:
- **Use online** — no installation, hosted limit applies.
- **Run locally** — clone the repository on macOS/Linux/Windows.
- **Windows + Joint-RPCA** — use WSL2 or Docker.
- **Docker** — optional reproducible local environment.

Link users to the README for full installation details rather than putting long shell setup procedures directly on Home.

### 17.6 README organization recommendation

Prefer an installation hierarchy such as:

```text
Using CorrOmics
├── Option 1 — Web app
├── Option 2 — Run locally from Git clone
│   ├── macOS
│   ├── Linux
│   └── Windows
├── Option 3 — Windows WSL2 (Joint-RPCA)
└── Option 4 — Docker
```

Within "Run locally from Git clone", keep the clone command common and only branch where OS-specific commands genuinely differ.

Do not duplicate full installation instructions in multiple README sections.


### 17.7 Git prerequisite and installation help

Under the Git-clone route, add a short prerequisite note:

> **Prerequisite:** Git must be installed before running `git clone`.

Keep this secondary to the main workflow so users who already have Git are not slowed down.

Provide concise OS-specific instructions:

**macOS**
Check:
```bash
git --version
```

If Git is not installed, offer one of:
```bash
xcode-select --install
```

or, if Homebrew is already installed:
```bash
brew install git
```

Do not require Homebrew just to install Git.

**Ubuntu/Debian Linux**
Check:
```bash
git --version
```

If Git is not installed:
```bash
sudo apt update
sudo apt install git
```

**Windows**
Check in PowerShell or Git Bash:
```powershell
git --version
```

If Git is not installed, direct the user to the official **Git for Windows** installer and keep the wording simple. Do not require users to understand Git Bash unless it is needed for a later command.

After installation, verify:
```bash
git --version
```

Then proceed with:
```bash
git clone https://github.com/abzer005/Corromics.git
cd Corromics
```

Keep these instructions easy to find under a small "Git not installed?" subsection or equivalent progressive-disclosure element.


### 17.8 Windows downloadable executable

Do not remove or de-emphasize the existing downloadable Windows executable option already documented in the repository.

Treat it as a separate access route alongside Web app, Git clone, WSL2, and Docker.

Recommended user-facing framing:

**Windows executable — easiest local Windows option**
- No manual Python/Conda environment setup should be required for normal use.
- All currently supported native-Windows CorrOmics functions should remain available.
- Joint-RPCA is **disabled in the native Windows executable** unless the application is explicitly redesigned to delegate Joint-RPCA to a separate WSL2/Docker Gemelli worker.
- Do not imply that installing Gemelli inside WSL automatically makes it available to the native Windows executable. WSL is a separate Linux environment.
- If future support is added, it would require an explicit bridge from the Windows app to a WSL/Docker Gemelli worker; treat that as a separate future enhancement, not part of this restructuring pass.

Preserve:
- executable-related build/configuration files;
- installer/spec files;
- README sections describing executable download/use;
- GitHub workflow/build files required to produce the executable.

Do not classify executable-related files as redundant merely because the app can also be run through Git, WSL, or Docker.

Add the executable to the access-route overview and method-availability table. The table should make Joint-RPCA unavailable for the native Windows executable unless verified otherwise.

Suggested installation/access hierarchy:

```text
Using CorrOmics
├── Option 1 — Web app
├── Option 2 — Windows executable
├── Option 3 — Run locally from Git clone
│   ├── macOS
│   ├── Linux
│   └── Windows
├── Option 4 — Windows WSL2 (Joint-RPCA)
└── Option 5 — Docker
```

Keep the executable route especially simple for non-technical Windows users.
