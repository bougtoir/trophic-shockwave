PY ?= python3

.PHONY: all data analysis figures report audit simulate clean help

help:
	@echo "make data      - download raw data + build processed panel"
	@echo "make analysis  - run event studies, distributed-lag, negative controls"
	@echo "make simulate  - run synthetic-data validation"
	@echo "make figures   - generate all figures and tables"
	@echo "make report    - render DATA_AUDIT/GO_NO_GO summaries"
	@echo "make all       - full pipeline"

all: data analysis simulate figures report

data:
	$(PY) -m src.download.dryad_download
	$(PY) -m src.download.modis_extract
	$(PY) -m src.download.climate_download
	$(PY) -m src.download.movebank_check
	$(PY) -m src.preprocess.build_panel
	$(PY) -m src.preprocess.vegetation_anomalies
	$(PY) -m src.preprocess.shock_events
	$(PY) -m src.preprocess.connectivity

analysis:
	$(PY) -m src.models.estimation
	$(PY) -m src.models.sensitivity
	$(PY) -m src.models.irf
	$(PY) -m src.models.negative_controls
	$(PY) -m src.models.autocorr

simulate:
	$(PY) -m src.simulation.simulate

figures:
	$(PY) -m src.figures.make_figures

report:
	$(PY) -m src.utils.render_reports
