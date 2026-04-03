# Project Layout

- `GEopf/` and `HGEopf/`: optimization code.
- `matpower7.1/`: third-party MATPOWER dependency.
- `wind data/`: bulk wind input files used by the workflow.
- `data/inputs/`: structured project inputs such as spreadsheets, geospatial files, weather samples, and network files.
- `data/optimization_snapshots/`: stored `.mat` snapshots from optimization experiments.
- `data/checkpoints/`: runtime checkpoint files such as `stop1.mat` through `stop6.mat`.
- `figs/`: generated figures exported by plotting scripts.
- `docs/reference/`: reference PDFs and supporting documents.
- `docs/archive/`: older exported documents kept for record.
- `logs/`: text logs.
- `archive/autosave/`: MATLAB autosave `.asv` files moved out of the root.
