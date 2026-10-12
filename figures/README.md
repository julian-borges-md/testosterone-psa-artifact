# Figures

Every figure in the manuscript, its generating code, and the exact values plotted.

| File | Figure | Panels | Plotted values |
|---|---|---|---|
| `Figure1_conceptual_model` | Causal model: testosterone, androgen receptor signaling, PSA, biopsy, diagnosis | Schematic | `data/Figure1_conceptual_model_edges.tsv` (each edge and the finding that supports it) |
| `Figure2_saturation_test` | Saturation and free hormone tests | Forest plot, 8 rows | `data/Figure2_saturation_test_data.tsv` |
| `Figure3_detection_signature` | Three exposures against PSA, PRACTICAL and FinnGen | 3 panels | `data/Figure3_detection_signature_data.tsv` |
| `Figure4_trials_and_synthesis` | A: randomized trials. B: posterior by assumed dose | 2 panels | `data/Figure4A_trial_inputs.tsv`, `data/Figure4B_posterior_curves.tsv` |

Each figure is supplied as 300 dpi PNG and vector PDF.

## Regenerate

```bash
make figures          # redraws all four figures from results/
python3 figures/export_figure_data.py   # rewrites data/ from results/
```

Both read only from `results/`; no value is retyped. Figure 4B reads the pooled prior synthesis (`results/stage4/synthesis_pooled.tsv`), so the posterior curves match the manuscript's primary analysis rather than the PRACTICAL only sensitivity analysis.

## Design

matplotlib, DejaVu Sans, 300 dpi, no top or right spines. Colors are fixed per exposure across figures: bioavailable testosterone dark blue, total testosterone grey, SHBG orange, BMI green, IGF-1 purple. Detection paths are red, tumour biology paths dashed green.
