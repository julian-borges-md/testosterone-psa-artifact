#!/usr/bin/env bash
# Build the Word versions of the venue neutral documents from their markdown
# sources. Pandoc resolves the figure references against figures/, so the two
# display figures are embedded from the same PNG files the figure scripts write.
#
#   bash manuscript/build.sh [output directory]
#
# The default output directory is build/. Every document is reproducible from
# the committed source: a rebuild of an unmodified source yields a document
# whose extracted text is identical, figure for figure, to the one before it.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/build}"
mkdir -p "$OUT"
command -v pandoc >/dev/null || { echo "pandoc is required" >&2; exit 1; }
cd "$ROOT"
for name in manuscript_v6.0 supplementary_v6.0 STROBE-MR_checklist figure_alt_text; do
  pandoc "manuscript/$name.md" --resource-path=figures -o "$OUT/$name.docx"
  echo "$OUT/$name.docx"
done
