# Expanded Selection-Function Status

Release: v1.0.6-review

The expanded parent-complete catalogue has been evaluated with the Castro-Ginard et al. (2023) Gaia DR3 RVS selection function through `gaiaunlimited.selectionfunctions.DR3RVSSelectionFunction` in a WSL Python environment. The per-star table is `expanded_selection_function.fits`; it includes `sf_parent_count` and `sf_prior_dominated_n_lt10` columns for reproducing the low-parent-count audit. Summary statistics are in `expanded_selection_function_summary.csv` and `expanded_selection_function_summary.json`; the table-level `n<10` reproduction audit is in `expanded_selection_function_lown_audit.csv`.

These inverse-selection weights are retained as contextual observability diagnostics, not as a volume-complete Milky Way deprojection or headline population count.
