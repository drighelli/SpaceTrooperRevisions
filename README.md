# SpaceTrooperRevisions

Collection of R scripts for running experiments and transferring/evaluating
query-score (QS) coefficients between spatial transcriptomics datasets and
platforms (e.g., CosMx, Xenium). The repository contains utilities for data
preparation, coefficient transfer, evaluation, and ablation studies.

Usage

- Run individual scripts from an R session or with `Rscript`, for example:

    Rscript run_qs_general.R

Scripts

- `ablation_study.R`: Perform ablation experiments to measure component impact.
- `compare_qs_general.R`: Compare query-score results across methods or datasets.
- `create_cosmx_rds.R`: Convert CosMx output into an RDS object for analysis.
- `evaluate_qs_split_kfold.R`: Evaluate QS performance using k-fold cross-validation.
- `run_qs_general.R`: General pipeline to compute and run query-score analyses.
- `transfer_coeffs_cosmx_cosmx.R`: Transfer coefficients within CosMx experiments.
- `transfer_coeffs_cosmx_xenium.R`: Transfer coefficients from CosMx to Xenium.
- `transfer_coeffs_xenium_xenium.R`: Transfer coefficients within Xenium experiments.
- `transfer_qs_coefficients.R`: Utilities for transferring QS coefficients between datasets.

Notes

- See individual scripts for usage details, required inputs, and expected
	outputs.

License

- See repository metadata or contact the author for licensing details.
