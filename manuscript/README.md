# Manuscript

| File | Contents |
|---|---|
| `manuscript_v1.0.md` | Full manuscript source. Tables are `TABLE:` blocks, equations `EQ:` lines, figures `FIG:` lines, references `REF:` lines |
| `build_docx.js` | Builds the submission Word file with the four figures embedded |

```bash
cd manuscript && npm install docx && node build_docx.js
```

Every number in the manuscript is asserted by `make verify` at the repository root. Items marked "to be confirmed" need the author's input before submission: funding, competing interests, the AI use statement, and the Zenodo DOI.
