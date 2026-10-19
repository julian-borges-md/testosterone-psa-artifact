# RLO_TPC | Is the prostate cancer risk of testosterone a PSA artifact?
#
#   make verify        re-derive every asserted number from committed results
#   make matrix        re-derive every cell of Tables 2, 3 and 4 from SNP level data
#   make manuscript    check the manuscript against the committed matrix and run the
#                      claim consistency, count qualification and supplementary gates
#   make submission    the above plus the unresolved placeholder gate
#   make crossvalidate compare every estimator against a published package (needs R)
#   make selftest      prove the assertion modules still fail on perturbed input
#   make figures       regenerate the display figures from committed results
#   make docx          build the Word versions of the documents (needs pandoc)
#   make archive-snapshot TAG=v6.3 DOI=10.5281/zenodo.N
#                      record what the cited deposit contains, after cutting a release
#   make clean         remove generated artefacts
#
# verify and selftest need no network and no raw data. Both run in CI.
# Full reanalysis from raw summary statistics: docs/RUNBOOK.md

PY := python3

.PHONY: install verify matrix matrix-full manuscript submission crossvalidate selftest figures docx archive-snapshot authorship hooks clean

install:
	$(PY) -m pip install -r requirements.txt

MS := manuscript/manuscript_v6.0.md

verify:
	$(PY) src/verify.py
	$(PY) src/verify_matrix.py
	$(PY) src/check_manuscript.py $(MS)

manuscript:
	$(PY) src/check_manuscript.py $(MS)

submission:
	$(PY) src/check_manuscript.py $(MS) --submission

crossvalidate:
	Rscript src/stage8_crossvalidation/crossvalidate.R
	Rscript src/stage8_crossvalidation/raps_reference.R
	Rscript src/stage8_crossvalidation/raps_loss_class.R

matrix:
	$(PY) src/verify_matrix.py

matrix-full:
	$(PY) src/verify_matrix.py --full

selftest:
	$(PY) src/selftest.py

figures:
	$(PY) figures/make_figures.py
	$(PY) figures/make_figure_saturation_full.py
	$(PY) figures/export_figure_data.py

docx:
	bash manuscript/build.sh

archive-snapshot:
	$(PY) src/write_archive_snapshot.py $(TAG) $(DOI)

authorship:
	$(PY) src/check_authorship.py

hooks:
	bash tools/install_hooks.sh

clean:
	find . -name __pycache__ -type d -exec rm -rf {} +
