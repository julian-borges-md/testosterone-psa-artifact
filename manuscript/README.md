# Manuscript

| File | Contents |
|---|---|
| `manuscript_v6.0.md` | Main text, tables and figure legends. Markdown: pipe tables, `![](figure.png)` image references resolved against `figures/` |
| `supplementary_v6.0.md` | Supplementary methods, tables S1 to S10 and the two supplementary figures |
| `STROBE-MR_checklist.md` | Completed reporting checklist |
| `title_page.md` | Title, author, correspondence, word and display item counts |
| `cover_letter.md` | Cover letter, venue neutral |
| `figure_alt_text.md` | The alt text and legends for both display figures, lifted from the manuscript for a submission system that asks for them in a field. Not a second source: `src/check_manuscript.py` holds it word for word identical to the manuscript |
| `build.sh` | Builds the Word versions with pandoc |

```bash
bash manuscript/build.sh            # writes build/*.docx
```

Every number in these documents is asserted by `make verify` at the repository
root, and `make submission` additionally asserts that no unresolved placeholder
remains. The documents name the archival deposit the repository corresponds to,
and `make manuscript` fails if that citation does not resolve to a version tag.

The sources are venue neutral by design. Correspondence addressed to a specific
journal is kept with the submission package rather than here.
