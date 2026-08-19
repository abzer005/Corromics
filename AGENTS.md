# AGENTS.md — Local Codex Instructions for CorrOmics

## Scope
Work only inside the CorrOmics repository currently opened in Codex.

- Do not read, search, inspect, or modify files outside the repository root.
- Do not access unrelated folders in the user's home directory.
- Do not change the scientific purpose or major workflow of CorrOmics.
- Before making any major conceptual, architectural, scientific-method, deployment-model, or data-model change, explain the proposed change and ask for permission.
- Small bug fixes, wording fixes, UI cleanup, and implementation changes explicitly required by the active spec may proceed without a second confirmation.

## Primary task
Follow `CORROMICS_RESTRUCTURE_SPEC.md` as the implementation specification for this restructuring pass.

Before editing:
1. Read `README.md`.
2. Inspect the complete repository tree.
3. Compare the repository layout with the Wang Bioinformatics Lab `Template_Streamlit_App` structure already referenced in the spec.
4. Identify files that appear redundant, obsolete, duplicated, generated, platform-specific, or unclear.
5. Do not delete or relocate uncertain files automatically. Report them first.

## Preserve
Do not break:
- Hosted GNPS deployment.
- Docker deployment.
- macOS local Git-clone workflow.
- Linux local Git-clone workflow.
- Windows WSL2 workflow.
- Native Windows standard workflow where currently supported.
- Joint-RPCA/Gemelli worker support.
- Existing example datasets.
- Existing scientific methods unless the spec explicitly changes their configuration.
- CSV and GraphML exports.

## Scientific-method guardrails

### SparCC
For CorrOmics cross-omics metabolite–microbiome analysis:
- Use the existing vendored SparCC implementation.
- Use deterministic normalization rather than the Dirichlet count model:
  - `norm="normalize"`
  - outer `iter=1`
- Keep SparCC's internal exclusion iterations (`xiter`, currently/default 10) unless there is a documented bug.
- Remove the user-facing outer SparCC iteration selector.
- Do not claim that SparCC "performs CLR internally."
- Explain that SparCC uses log-ratio variance inference and is exploratory for heterogeneous metabolite–microbiome cross-omics data.
- Do not add `+1` to the already-normalized combined target dataframe inside SparCC.
- The second-omics preprocessing should provide a zero-safe relative-abundance option by adding a pseudocount before relative-abundance normalization when the user needs SparCC-compatible data.
- Do not silently change the meaning of the existing zero-preserving relative-abundance option.

### Joint-RPCA
- Do not change the scientific Joint-RPCA calculations unless a genuine implementation bug is demonstrated.
- Support a separate Gemelli worker interpreter via `GEMELLI_WORKER_PYTHON`.
- Do not require Gemelli to be installed in the same Python environment as Streamlit.
- Make dependency/error guidance OS-aware.
- Native Windows should recommend WSL2 or Docker for Joint-RPCA.
- macOS errors must not incorrectly say Linux/WSL is required.
- Linux/WSL should receive Linux-specific guidance.
- Docker should honor its configured Gemelli worker.

## Terminology
Use **association** for the general cross-omics workflow, result collections, and generic UI labels.
Use **correlation** when referring specifically to:
- Pearson correlation
- Spearman correlation
- SparCC correlation coefficients
- a mathematically defined correlation coefficient

Do not rename scientific method names.

## UI rules
- Preserve the current overall app concept and navigation.
- Prefer compact tabs where explicitly requested in the spec.
- Avoid adding large new warning blocks.
- Warnings should be actionable, method-specific, and concise.
- Do not duplicate the same warning in multiple places when one reusable helper can provide it.

## Deployment-mode rules
- Hosted/restricted GNPS mode: maximum 1,000,000 cross-omics pairs.
- Local mode: maximum 10,000,000 cross-omics pairs.
- Do not infer "local" merely because the app runs inside Docker.
- Prefer an explicit deployment/environment flag over unreliable runtime heuristics.
- Keep local Docker and hosted Docker distinguishable.

## Repository cleanup
When auditing files:
- Categorize each as runtime, deployment, development/build, documentation, example data, third-party/vendor, generated artifact, or candidate for removal.
- Do not delete candidate files until the user approves the removal list.
- Preserve third-party license/copyright notices.
- Do not rewrite vendored third-party code unless required for compatibility and clearly documented.

## Licensing
Do not make a final license change in this pass unless explicitly approved.
- The intended top-level CorrOmics license is currently BSD 3-Clause, pending final review.
- Vendored third-party components may have different licenses and must retain their applicable notices.
- Produce a licensing audit before removing or replacing any license file.

## Versioning
The public repository currently has no GitHub Releases/tags and `pyproject.toml` does not define an application version.
For this restructuring:
- Introduce a simple, single-source application version if practical.
- Use Semantic Versioning (`MAJOR.MINOR.PATCH`) unless existing code provides a better established convention.
- Do not create/push a GitHub release or Git tag automatically.
- Prepare the repository so the user can make this restructuring the first formal release.
- Add/update a changelog or release-notes mechanism if absent.
- Report the proposed version and rationale before tagging.

## Validation
After changes:
- Run syntax/import checks.
- Run available tests.
- Launch/test the main Streamlit app where feasible.
- Check key pages affected by the change.
- Verify local vs hosted limit logic.
- Verify SparCC output on a small positive test matrix.
- Verify Joint-RPCA worker resolution logic without changing scientific calculations.
- Verify README commands for macOS/Linux/WSL/Docker are internally consistent.

## Final report
At completion, provide:
1. Files changed.
2. Files proposed for removal, with reasons.
3. Files intentionally kept despite appearing redundant, with reasons.
4. Scientific behavior changed.
5. Scientific behavior intentionally unchanged.
6. Deployment behavior changed.
7. Documentation changes.
8. Tests/checks run and results.
9. Remaining issues/limitations.
10. Proposed release version.
11. `git diff --stat`.
12. Important diffs requiring manual review.

Do not commit, tag, push, publish a release, or delete uncertain files without explicit user approval.


## Installation UX
Treat installation documentation as scientist-facing documentation, not developer-facing documentation.

- Present the Web app as the easiest no-install route.
- Present Git clone as the main local route and recognize that Git clone works on macOS, Linux, and Windows when Git is installed.
- Present WSL2 specifically as a Windows compatibility route, especially for Joint-RPCA.
- Present Docker as an optional cross-platform reproducible/isolation route, not as a Windows-only feature.
- Keep advanced Conda/Docker/environment details behind the simple choice guidance.
- Verify method availability per access route before documenting it.
- Do not expose implementation jargon unless it is necessary for troubleshooting.


- Under Git-clone instructions, include a short "Git not installed?" path for macOS, Ubuntu/Debian Linux, and Windows, plus `git --version` verification. Keep it secondary to the main install flow.


## Windows executable
- Preserve the existing Windows executable/installer route and all files needed to build/distribute it.
- Do not mark executable-related files as redundant without proving they are unused.
- Treat the executable as the simplest local Windows option for non-technical users.
- Joint-RPCA remains unavailable in the native Windows executable unless there is an explicit implementation that delegates to a WSL2/Docker Gemelli worker.
- Do not claim that installing Gemelli inside WSL automatically enables Joint-RPCA in the native Windows executable.
