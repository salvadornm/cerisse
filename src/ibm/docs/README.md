# IBM Documentation

This directory stores maintained documentation for the immersed boundary module in `src/ibm`.

## Primary Documents

- **`ibm_complete.md`** — Single-file Markdown reference covering: code layout, CGAL→BVH type mapping, initialization flow (Stages A–D), and runtime path.
- **`latex/ibm_complete.tex`** — LaTeX version of the same consolidated reference. Compile with `pdflatex ibm_complete.tex`.

## Legacy (individual section files)

The following files are kept as-is for reference but `ibm_complete.md/tex` is the canonical source:

- `latex/ibm_source_reference.tex` — file layout table only
- `latex/ibm_initialization.tex` — initialization stages only
- `latex/ibm_bvh_type_mapping.tex` — CGAL→BVH type mapping table only
- `markdown/ibm_bvh_type_mapping.md` — Markdown type mapping table only

Generated artifacts (PDF, aux, log, etc.) are intentionally not stored here.
