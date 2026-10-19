"""RLO_TPC | Manuscript concordance and claim consistency check

Parses the manuscript markdown and asserts, on every push, that the document
agrees with the computation and with itself.

1. Concordance. Every cell printed in Tables 2, 3 and 4, and every numeric claim
   in the abstract that corresponds to a table cell, matches
   results/stage7_matrix/estimator_matrix.tsv. An abstract and a table can
   disagree silently whenever the matrix is recomputed and only one of them is
   updated, so neither is trusted against the other; both are checked against
   the committed source.

2. Claim consistency. No sentence in the abstract, the Key Messages or the
   Discussion opening asserts a falsification that the limitations table denies.
   Compressing a manuscript strips qualifications before it strips claims,
   because a hedge costs words and an assertion does not, so this is run as a
   step rather than left to a careful reading.

3. Count qualification. Any count that rests on the MR RAPS row states the loss
   setting it depends on, and Supplementary Table S6 is cited from the main
   text. Which counts require it is computed from the committed loss
   sensitivity tables rather than listed here, so a regeneration that makes a
   stable count setting dependent begins failing this gate unprompted.

4. Supplementary citation. Every supplementary table is cited at least once
   from the main text. An uncited supplementary table is in practice an unread
   one.

5. Figure concordance. Every plotted value equals its committed matrix cell,
   and no figure labels an estimator with a method qualifier the manuscript
   does not use.

6. Companion documents. The reporting checklist, the cover letter, CITATION.cff
   and the README each carry the manuscript's current title and none of the
   claims the manuscript has retracted.

7. Venue neutrality. The tracked tree names no submission target. Tokens are in
   docs/venue_denylist.txt rather than here, so that this module is scanned too.

8. Archive citation. The tag the deposit corresponds to is cited and has the
   shape of a version tag, so a prose edit cannot leave the citation
   grammatical and pointing at nothing.

9. Falsification polarity. In the manuscript and in every document that travels
   with it, a falsification token near the saturation subject carries a negation
   before it. Gate 6 lists retracted phrases and is defeated by paraphrase; this
   states the property, which is what the limitations table actually asserts.

10. Absence limit. Every document that summarises the study carries the sentence
   that the evidence does not establish absence of effect. A summary that keeps
   the detection half of the conclusion and drops the limit is a stronger paper
   than the one submitted.

11. Archive content. The documents the cited deposit is supposed to contain are
   unchanged since it was cut, recorded in docs/archive_snapshot.tsv. Gate 8
   checks that the citation has the shape of a tag; this checks that the tag
   holds this manuscript, which is the failure that actually occurred.

Journal limits are asserted on a documented basis. With --submission, the module
additionally asserts that no unresolved placeholder remains; that gate is a
submission precondition and deliberately does not run on every push, so that a
known pending item does not mask a real regression in continuous integration.

usage: python3 src/check_manuscript.py path/to/manuscript.md [--submission]
"""

import hashlib
import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = sys.argv[1]
SUBMISSION = "--submission" in sys.argv[2:]
M = pd.read_csv(os.path.join(ROOT, "results", "stage7_matrix", "estimator_matrix.tsv"), sep="\t")
C = pd.read_csv(os.path.join(ROOT, "results", "stage7_matrix",
                             "contrast_biot_minus_tt.tsv"), sep="\t")
text = open(MS).read()

# Journal limits are a property of the venue, not of the manuscript, so they are
# committed as data rather than written into this module. A null means the
# journal publishes no limit for that item, and the check is then not run.
LIMITS = {"main_text_words": None, "abstract_words": 250, "title_characters": None}
_lim = os.path.join(ROOT, "docs", "journal_limits.tsv")
if os.path.exists(_lim):
    with open(_lim) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                k, v = line.split("\t")[:2]
                LIMITS[k.strip()] = int(v) if v.strip().isdigit() else None

ROWMAP = {
    "IVW": "IVW",
    "MR Egger": "MR Egger",
    "Weighted median": "Weighted median",
    "Weighted mode": "Weighted mode",
    "IVW, Q outliers removed": "IVW, Q outliers removed",
    "MR cML BIC, data perturbation SE": "MR cML BIC (DP SE)",
    "MR RAPS Huber (mr.raps 0.2)": "MR RAPS Huber",
}
T2COLS = ["biot->psa", "tt->psa", "shbg->psa"]
T4COLS = ["biot->prac", "tt->prac", "shbg->prac", "biot->fg"]
CELL = re.compile(r"\*{0,2}(-?\d+\.\d+)\*{0,2}\s*\((-?\d+\.\d+) to (-?\d+\.\d+)\)")

fails, n = [], 0


def table_rows(header_fragment):
    out = []
    grab = False
    for line in text.splitlines():
        if header_fragment in line:
            grab = True
            continue
        if grab:
            if line.startswith("|"):
                cells = [c.strip() for c in line.strip("|").split("|")]
                out.append(cells)
            elif out:
                break
    return out


def cmp_cell(tag, printed, beta, se, scale):
    global n
    n += 1
    lo, hi = beta - 1.96 * se, beta + 1.96 * se
    if scale == "or":
        beta, lo, hi = np.exp(beta), np.exp(lo), np.exp(hi)
    m = CELL.search(printed)
    if not m:
        fails.append(f"{tag}: could not parse printed cell {printed!r}")
        return
    pb, plo, phi = (float(x) for x in m.groups())
    for lab, got, exp in (("point", pb, beta), ("lower", plo, lo), ("upper", phi, hi)):
        if abs(got - exp) > 5.1e-4:
            fails.append(f"{tag} {lab}: manuscript {got} vs committed {exp:.4f}")
    bold = printed.strip().startswith("**")
    sig = (lo > 0 and hi > 0) or (lo < 0 and hi < 0) if scale == "log" else (lo > 1) or (hi < 1)
    if bold != bool(sig):
        fails.append(f"{tag}: bold marking {bold} does not match interval excluding the null {sig}")


for frag, cols, scale in [("| Estimator | Bioavailable testosterone |", T2COLS, "log"),
                          ("| Estimator | Bioavailable T, PRACTICAL |", T4COLS, "or")]:
    for cells in table_rows(frag):
        row = cells[0].strip()
        if row.startswith("---") or row not in ROWMAP:
            continue
        est = ROWMAP[row]
        for i, a in enumerate(cols):
            rec = M[(M.analysis == a) & (M.estimator == est)].iloc[0]
            cmp_cell(f"{'T2' if scale == 'log' else 'T4'} {a} {est}", cells[i + 1],
                     rec.beta, rec.se, scale)

for cells in table_rows("| Estimator | Bioavailable | Total | Difference (SE) |"):
    row = cells[0].strip()
    if row.startswith("---") or row not in ROWMAP:
        continue
    rec = C[C.estimator == ROWMAP[row]].iloc[0]
    m = re.match(r"([+-]?\d+\.\d+)\s*\(([\d.]+)\)", cells[3].strip())
    if not m:
        fails.append(f"T3 {row}: could not parse difference cell {cells[3]!r}")
        continue
    n += 1
    if abs(float(m.group(1)) - rec.difference) > 5.1e-5:
        fails.append(f"T3 {row} difference: manuscript {m.group(1)} vs committed {rec.difference:.4f}")
    n += 1
    if abs(float(m.group(2)) - rec.se) > 5.1e-5:
        fails.append(f"T3 {row} se: manuscript {m.group(2)} vs committed {rec.se:.4f}")

# ---- abstract against the matrix ----
abstract = text.split("**Results.**")[1].split("**Discussion.**")[0]


def abstract_check(label, pattern, analysis, estimator, scale):
    global n
    n += 1
    m = re.search(pattern, abstract)
    if not m:
        fails.append(f"abstract: pattern for {label} not found")
        return
    rec = M[(M.analysis == analysis) & (M.estimator == estimator)].iloc[0]
    exp = np.exp(rec.beta) if scale == "or" else rec.beta
    printed = m.group(1)
    got = float(printed)
    # tolerate rounding at the precision the abstract prints, and no more
    dp = len(printed.split(".")[1])
    if abs(got - exp) > 0.5 * 10 ** (-dp) + 1e-9:
        fails.append(f"abstract {label}: {got} vs committed {exp:.4f} "
                     f"(printed to {dp} decimal places)")


abstract_check("IVW on PSA", r"inverse variance weighted (\d+\.\d+) natural log units",
               "biot->psa", "IVW", "log")
abstract_check("PRACTICAL OR", r"PRACTICAL \(OR (\d+\.\d+),", "biot->prac", "IVW", "or")
abstract_check("FinnGen OR", r"FinnGen \(OR (\d+\.\d+),", "biot->fg", "IVW", "or")

# counts of estimators whose interval excludes the null, as claimed in the abstract
n += 1
k = int(M[(M.analysis == "biot->psa") & M.excludes_null].shape[0])
if "six of seven estimators" not in abstract or k != 6:
    fails.append(f"abstract: bioavailable testosterone significant in {k} of 7")

# counts asserted against the Results text rather than the abstract, which the
# word limit no longer leaves room for
results = text.split("**Results**")[1].split("**Discussion**")[0]
for label, pattern, got in [
    ("total testosterone on PSA", "significantly positive on PSA under two of seven",
     int(M[(M.analysis == "tt->psa") & M.excludes_null & (M.beta > 0)].shape[0])),
    ("SHBG on PSA", "SHBG lowers PSA in six of seven",
     int(M[(M.analysis == "shbg->psa") & M.excludes_null & (M.beta < 0)].shape[0])),
    ("bioavailable T on PRACTICAL", "in PRACTICAL in five of seven",
     int(M[(M.analysis == "biot->prac") & M.excludes_null].shape[0])),
    ("total testosterone on PRACTICAL", "significantly positive on PRACTICAL in two of seven",
     int(M[(M.analysis == "tt->prac") & M.excludes_null & (M.beta > 0)].shape[0])),
]:
    n += 1
    want = {"total testosterone on PSA": 2, "SHBG on PSA": 6,
            "bioavailable T on PRACTICAL": 5, "total testosterone on PRACTICAL": 2}[label]
    if pattern not in results or got != want:
        fails.append(f"results: {label} is {got} of 7, text pattern {pattern!r} "
                     f"{'present' if pattern in results else 'ABSENT'}")
n += 1
sig = int((C.p < 0.05).sum())
if sig != 6 or "intervals excluding zero under six" not in results:
    fails.append(f"results: contrast significant under {sig} of 7 estimators")
n += 1
pmax = C[C.p < 0.05].p.max()
if abs(pmax - 0.0047) > 5e-5 or "to 0.0047" not in results:
    fails.append(f"results: largest significant contrast p is {pmax:.4f}")

# ---------------------------------------------------------------------------
# claim consistency, count qualification, supplementary citation, limits
# ---------------------------------------------------------------------------


def region(start, end):
    i = text.index(start)
    j = text.index(end, i + len(start))
    return text[i:j]


abstract_region = region("**Abstract**", "**Keywords:**")
# The Key Messages element belongs to a journal that requires it. Where the
# target journal has no such element the region is empty, and the claim
# consistency gate runs over the regions that do exist rather than failing.
keymsg_region = (region("**Key Messages**", "**Introduction**")
                 if "**Key Messages**" in text else "")
disc_region = region("**Discussion**", "The dogma holds")
results_region = region("**Results**", "**Discussion**")
main_text = region("**Introduction**",
                   "**Acknowledgments**" if "**Acknowledgments**" in text
                   else "**Declarations**")

# -- gate 2: claim consistency ---------------------------------------------
# The limitations table denies that a linear instrumental variable average can
# falsify a threshold model, so no prominent region may assert that it did.
# Checked mechanically because a length cut removes qualifications before it
# removes claims.
DENIED = [
    ("falsif", "asserts a falsification"),
    ("is not saturated", "asserts the prostate is not saturated"),
    ("prostate is not saturated", "asserts the prostate is not saturated"),
]
for label, reg in (("abstract", abstract_region), ("Key Messages", keymsg_region),
                   ("Discussion opening", disc_region)):
    for needle, what in DENIED:
        n += 1
        if needle in reg.lower():
            fails.append(f"claim consistency: {label} {what} "
                         f"(matched {needle!r}), which the limitations row "
                         f"'A linear IV average cannot falsify a threshold model' denies")

# -- gate 3: count qualification --------------------------------------------
# Two properties are needed and a phrase list gives neither.
#
# Which statements need qualification is computed from the committed loss
# sensitivity tables: the MR RAPS row is substituted under each loss and a
# statement needs qualification exactly when its count moves. Nothing here names
# a particular exposure, so a regeneration that makes a stable count setting
# dependent starts failing this gate without anyone remembering to extend it.
#
# And the qualification must be LOCAL to the statement. A region wide word
# search passes by accident, because a count word that qualifies one statement
# sits in the same region as another it does not qualify. Each computed
# statement therefore declares the phrase that prints it, the alternative must
# appear within WINDOW characters of that phrase together with a loss name, and
# a statement whose phrase cannot be found fails rather than passing silently.
QUALIFIERS = ("l2", "huber", "tukey")
WINDOW = 300
WORD = {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four",
        5: "five", 6: "six", 7: "seven"}
LOSSDIR = os.path.join(ROOT, "results", "stage8_crossvalidation")
LOSS_E = os.path.join(LOSSDIR, "raps_loss_by_exposure.tsv")
LOSS_C = os.path.join(LOSSDIR, "raps_loss_contrast.tsv")

# each entry: key, the phrase in the Results that prints the statement
ANCHOR = {
    ("biot->psa", "pos"): "excluding the null in six of seven estimators",
    ("tt->psa", "pos"): "significantly positive on PSA under two of seven",
    ("tt->psa", "neg"): "significantly negative under MR Egger",
    ("shbg->psa", "neg"): "SHBG lowers PSA in six of seven",
    ("contrast", "sig"): "intervals excluding zero under six",
    ("contrast", "pmax"): "to 0.0047",
}


def local_window(anchor):
    i = results_region.find(anchor)
    if i < 0:
        return None
    return results_region[i:i + WINDOW].lower()


def require_local(key, adopted, got):
    """Require every alternative count to be stated beside its own statement."""
    global n
    alternatives = sorted({v for v in got.values() if v != adopted})
    n += 1
    anchor = ANCHOR[key]
    win = local_window(anchor)
    if win is None:
        fails.append(f"count qualification: cannot locate the statement for {key} in the "
                     f"Results (expected the phrase {anchor!r}), so the gate cannot "
                     f"confirm whether it is qualified")
        return
    if not alternatives:
        return
    missing = [WORD[v] for v in alternatives if WORD[v] not in win]
    if missing or not any(q in win for q in QUALIFIERS):
        fails.append(f"count qualification: {key} is {adopted} under the adopted MR RAPS "
                     f"setting and {alternatives} under the others, and the statement at "
                     f"{anchor!r} does not state "
                     f"{', '.join(missing) or 'the loss setting'} within {WINDOW} characters")


n += 1
if not (os.path.exists(LOSS_E) and os.path.exists(LOSS_C)):
    fails.append("count qualification: the committed loss sensitivity tables are "
                 "absent, so the class level gate cannot run")
else:
    LE = pd.read_csv(LOSS_E, sep="\t")
    LCON = pd.read_csv(LOSS_C, sep="\t")
    losses = sorted(set(LE.loss))

    for analysis, direction in (("biot->psa", "pos"), ("tt->psa", "pos"),
                                ("tt->psa", "neg"), ("shbg->psa", "neg")):
        got = {}
        for lf in losses:
            r = LE[(LE.analysis == analysis) & (LE.loss == lf)]
            if len(r) != 1:
                continue
            r = r.iloc[0]
            sub = M[M.analysis == analysis].set_index("estimator").copy()
            sub.loc["MR RAPS Huber", ["beta", "lo", "hi"]] = [r.beta, r.lo, r.hi]
            got[lf] = int(((sub.lo > 0) & (sub.hi > 0)).sum() if direction == "pos"
                          else ((sub.lo < 0) & (sub.hi < 0)).sum())
        if got:
            require_local((analysis, direction), got.get("huber"), got)

    csig, cpmax = {}, {}
    for lf in losses:
        r = LCON[LCON.loss == lf]
        if len(r) != 1:
            continue
        r = r.iloc[0]
        cc = C.set_index("estimator").copy()
        cc.loc["MR RAPS Huber", ["difference", "se", "p"]] = [r.difference, r.se, r.p]
        sig = cc[cc.p < 0.05]
        csig[lf] = len(sig)
        cpmax[lf] = "%.4f" % sig.p.max()
    if csig:
        require_local(("contrast", "sig"), csig.get("huber"), csig)
    n += 1
    if len(set(cpmax.values())) > 1:
        win = local_window(ANCHOR[("contrast", "pmax")])
        if win is None or not any(q in win for q in QUALIFIERS):
            fails.append(f"count qualification: the printed contrast p upper bound differs "
                         f"across losses ({cpmax}) and the Results state only one value")

n += 1
if "Table S6" not in main_text:
    fails.append("count qualification: Supplementary Table S6, which governs the "
                 "counts that rest on the MR RAPS row, is not cited from the main text")

# -- gate 4: supplementary citation -------------------------------
for i in range(1, 11):
    n += 1
    if not re.search(r"\bS%d\b" % i, main_text):
        fails.append(f"supplementary citation: Table S{i} is never cited from the main text")

# -- gate 5: figure concordance ---------------------------------------------
# The figures are derived from the matrix and read by the same reviewer as the
# tables, and a module that parses only the manuscript does not cover them.
# Two things are asserted: plotted values equal their committed matrix cells,
# and any parenthetical method qualifier on a figure label is language the
# manuscript itself uses, so a figure cannot describe an estimator in terms the
# paper does not.
FIGDATA = os.path.join(ROOT, "figures", "data",
                       "Figure1_saturation_full_matrix_data.tsv")
PANEL = {"Bioavailable testosterone": "biot->psa",
         "Total testosterone": "tt->psa", "SHBG": "shbg->psa"}
FIGEST = {"IVW": "IVW", "MR Egger": "MR Egger", "Weighted median": "Weighted median",
          "Weighted mode": "Weighted mode",
          "IVW, Q outliers removed": "IVW, Q outliers removed",
          "MR cML (data perturbation SE)": "MR cML BIC (DP SE)",
          "MR RAPS (mr.raps 0.2)": "MR RAPS Huber"}

n += 1
if not os.path.exists(FIGDATA):
    fails.append("figure concordance: the exported figure data is absent, so the "
                 "plotted values cannot be checked against the committed matrix")
else:
    F = pd.read_csv(FIGDATA, sep="\t")
    for _, row in F.iterrows():
        n += 1
        a, e = PANEL.get(row.panel), FIGEST.get(row.estimator)
        if a is None or e is None:
            fails.append(f"figure concordance: the figure plots panel {row.panel!r} "
                         f"estimator {row.estimator!r}, which does not map to a "
                         f"committed matrix cell")
            continue
        rec = M[(M.analysis == a) & (M.estimator == e)].iloc[0]
        for lab, got, exp in (("beta", row.beta, rec.beta), ("lo", row.lo, rec.lo),
                              ("hi", row.hi, rec.hi)):
            if abs(got - exp) > 1e-9:
                fails.append(f"figure concordance: {row.panel} {row.estimator} {lab} "
                             f"plotted {got} against committed {exp}")
    # a parenthetical method qualifier on a figure label must be language the
    # manuscript itself uses for that estimator
    for label in sorted(set(F.estimator)):
        m = re.search(r"\(([^)]+)\)", label)
        if not m:
            continue
        n += 1
        if m.group(1) not in text:
            fails.append(f"figure concordance: the figure labels an estimator "
                         f"{label!r}, and the qualifier {m.group(1)!r} appears "
                         f"nowhere in the manuscript")

# -- gate 6: companion documents --------------------------------------------
# The reporting checklist, the cover letter, CITATION.cff and the README travel
# with the manuscript and are read by the same editor. Each must carry the
# current title, so a superseded description cannot travel with the submission,
# and none may assert anything the manuscript has retracted. The denylist holds
# assertions rather than the word \"falsify\", because the correct usage of that
# word is a denial: the limitations row and the checklist both use it that way,
# and a bare word match would flag the sentences that state the position right.
COMPANIONS = ["manuscript/STROBE-MR_checklist.md", "manuscript/cover_letter.md",
              "CITATION.cff", "README.md"]
RETRACTED = ["falsifying the saturation prediction",
             "falsifies the central prediction",
             "makes a falsifiable prediction",
             "prespecified falsification criterion",
             "a falsification criterion",
             "prostate is not saturated",
             "responds to bioavailable, not total"]
title_line = text.strip().splitlines()[0].replace("**", "").strip()

for rel in COMPANIONS:
    path = os.path.join(ROOT, rel)
    n += 1
    if not os.path.exists(path):
        fails.append(f"companion documents: {rel} is missing, and it travels with "
                     f"the submission")
        continue
    comp = open(path).read()
    n += 1
    if title_line not in comp:
        fails.append(f"companion documents: {rel} does not carry the manuscript's "
                     f"current title, so it describes a superseded version")
    for phrase in RETRACTED:
        n += 1
        if phrase in comp.lower() if phrase.islower() else phrase in comp:
            fails.append(f"companion documents: {rel} asserts {phrase!r}, which the "
                         f"manuscript has retracted")

# -- gate 9: falsification polarity -----------------------------------------
# Gate 6 holds a list of retracted phrases, which is an enumeration and is
# therefore defeated by paraphrase. It carried "prespecified falsification
# criterion" while the README said "Prespecified falsification rule met" and
# CITATION.cff said "falsifies the androgen saturation model", so an assertion
# the manuscript disclaims travelled in the public record and in the archival
# metadata. This gate states the property instead of listing its instances.
#
# A falsification token near the saturation subject is an assertion unless it
# is negated, and in every correct usage the negation PRECEDES the token: "does
# not claim to falsify", "cannot falsify". Only the preceding text is searched,
# so a trailing clause that happens to contain a negative word cannot excuse an
# affirmative one, which is how "The saturation falsification ... do not depend
# on it" would otherwise pass. Tokens with no saturation subject nearby are not
# this gate's business: the falsifiability of the verification harness is a
# different and correct use of the word.
#
# The search runs to the start of the sentence or table cell holding the token
# and no further. A fixed character window was tried first and leaked at once:
# the paragraph above the findings table ends "does not claim it", so every
# claim in the first row of that table inherited a negation from a different
# sentence. A negation only qualifies the clause it is in.
# Coverage is the documents a reader reaches, which is wider than the ones that
# travel with the submission: the repository front page, the manifest in docs
# and the figure scripts all carried the claim too, and the figure scripts
# print their text into a published image.
#
# prespecification/ is deliberately NOT covered. It records what was decided
# before estimation, in the words used then, including a decision rule written
# as "Saturation FALSIFIED if ...". Editing that to match what the analysis
# later supported would not be correcting a claim, it would be falsifying the
# record of the claim. A prespecification is evidence, not prose, and this gate
# must never be the reason someone rewrites one.
POLARITY = ["manuscript/manuscript_v6.0.md", "manuscript/supplementary_v6.0.md",
            "manuscript/STROBE-MR_checklist.md", "manuscript/cover_letter.md",
            "CITATION.cff", "README.md"]
for _d, _ext in (("docs", ".md"), ("figures", ".py")):
    _full = os.path.join(ROOT, _d)
    if os.path.isdir(_full):
        POLARITY += [os.path.join(_d, f) for f in sorted(os.listdir(_full))
                     if f.endswith(_ext)]
TOKEN = re.compile(r"falsif\w*", re.I)
SUBJECT = re.compile(r"saturation|saturated|threshold model|central prediction", re.I)
NEGATOR = re.compile(r"\bnot\b|\bno\b|\bcannot\b|\bnever\b|\bdenies\b|\bdenied\b|"
                     r"\brather than\b", re.I)
SUBJ, UNIT_MAX = 200, 300
BOUNDARY = re.compile(r"[|\n]|\. |: |; ")


def clause_before(doc, i):
    """The text from the start of the sentence or table cell up to i."""
    seg = doc[max(0, i - UNIT_MAX):i]
    last = 0
    for b in BOUNDARY.finditer(seg):
        last = b.end()
    return seg[last:]


for rel in POLARITY:
    path = os.path.join(ROOT, rel)
    if not os.path.exists(path):
        continue
    doc = open(path, encoding="utf-8", errors="ignore").read()
    for mt in TOKEN.finditer(doc):
        near = doc[max(0, mt.start() - SUBJ):mt.end() + SUBJ]
        if not SUBJECT.search(near):
            continue
        n += 1
        if not NEGATOR.search(clause_before(doc, mt.start())):
            fails.append(f"falsification polarity: {rel} asserts "
                         f"{doc[max(0, mt.start() - 60):mt.end() + 40].strip()!r} "
                         f"with no negation before the claim, and the manuscript "
                         f"holds that a linear IV average cannot refute a "
                         f"threshold model")

# -- gate 10: the absence limit travels with every summary ------------------
# The manuscript's conclusion is bounded on both sides: part of the apparent
# risk runs through detection, AND the evidence does not establish absence of
# effect. The first half was itself overstated as "much" until 19 October and is
# now held to a bound by the share claim gate below. A summary that keeps the
# first half and drops the second is a
# stronger paper than the one submitted, and the data availability statement
# sends readers to these documents. The positive requirement is the property
# doing the work here; the two denied strings are the historical forms this
# repository actually carried, kept because an exact form costs nothing to
# check once it is known.
LIMIT = "does not establish absence of effect"
CARRIES_LIMIT = ["manuscript/manuscript_v6.0.md", "CITATION.cff", "README.md"]
ABSENCE_FORMS = ["raises psa, not prostate cancer risk",
                 "no meaningful effect of testosterone on prostate cancer"]

for rel in CARRIES_LIMIT:
    path = os.path.join(ROOT, rel)
    n += 1
    if not os.path.exists(path):
        fails.append(f"absence limit: {rel} is missing")
        continue
    doc = open(path, encoding="utf-8", errors="ignore").read()
    if LIMIT not in doc:
        fails.append(f"absence limit: {rel} summarises the study without the "
                     f"sentence that the evidence {LIMIT!r}, which the abstract "
                     f"carries and which bounds the conclusion")
    for form in ABSENCE_FORMS:
        n += 1
        if form in doc.lower():
            fails.append(f"absence limit: {rel} asserts {form!r}, which states "
                         f"an absence of effect that the manuscript declines to "
                         f"claim")

# -- gate 7: the repository names no submission target ----------------------
# This repository holds the analysis and the manuscript and is venue neutral;
# submission correspondence, which names a journal by nature, is kept outside
# it. The scan covers the whole tracked tree rather than the manuscript alone,
# because a venue name reaches a figure docstring or a module comment as readily
# as the prose. The tokens live in docs/venue_denylist.txt rather than here, so
# that this module is itself scanned; the denylist is the one exempt file.
DENYLIST = os.path.join(ROOT, "docs", "venue_denylist.txt")
SKIP_DIRS = {".git", "__pycache__"}
SKIP_EXT = {".png", ".pdf", ".gz", ".tsv", ".csv", ".bundle"}

n += 1
if not os.path.exists(DENYLIST):
    fails.append("submission target: docs/venue_denylist.txt is absent, so the "
                 "venue neutrality of the repository cannot be checked")
else:
    literals, patterns = [], []
    for raw in open(DENYLIST):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if len(line) > 2 and line.startswith("/") and line.endswith("/"):
            patterns.append(re.compile(line[1:-1]))
        else:
            literals.append(line)
    hits = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fname in filenames:
            full = os.path.join(dirpath, fname)
            rel = os.path.relpath(full, ROOT)
            if rel == os.path.join("docs", "venue_denylist.txt"):
                continue
            if os.path.splitext(fname)[1].lower() in SKIP_EXT:
                continue
            try:
                body = open(full, encoding="utf-8", errors="ignore").read()
            except OSError:
                continue
            for token in literals:
                if token in body:
                    hits.append("%s names %r" % (rel, token))
            for pat in patterns:
                if pat.search(body):
                    hits.append("%s matches /%s/" % (rel, pat.pattern))
    if hits:
        fails.append("submission target: this repository is venue neutral, and "
                     + "; ".join(sorted(set(hits))))

# -- gate 8: the archive citation is well formed ----------------------------
# The data availability statement and the Methods cite a git tag, which is what
# tells a reader which state of the repository the archival deposit corresponds
# to. A prose edit can leave that citation grammatical and meaningless, so the
# shape is asserted rather than assumed: the token after "at tag" must look like
# a version tag and not like prose.
n += 1
tags = [m.rstrip(".,;:") for m in re.findall(r"(?<![A-Za-z])at tag (\S+)", text)]
if not tags:
    fails.append("archive citation: the manuscript does not cite the tag the "
                 "archival deposit corresponds to")
for tag in tags:
    n += 1
    if not re.fullmatch(r"v\d+\.\d+(\.\d+)?", tag):
        fails.append(f"archive citation: {tag!r} is cited as a tag but is not a "
                     f"version tag, so the citation points at nothing")

# -- gate 11: the cited archive contains this manuscript --------------------
# Gate 8 asserts that the tag citation has the shape of a tag. It does not
# assert that the tag contains this manuscript, and twice now it has not: a tag
# was cut, the manuscript was then corrected, and the citation went on naming a
# deposit holding a different document, the second time with a different title.
# A reader who follows the DOI is entitled to the document that cites it. The
# condition is recorded as data by src/write_archive_snapshot.py and checked
# here.
#
# Hashes are taken with every DOI and version tag masked. The commit that
# updates the citation after a tag is cut cannot be inside that tag, so it is
# the one permitted difference between the deposit and the submitted document,
# and masking is what permits exactly it and nothing else. The verification
# modules are deliberately outside the coverage: this gate ships in that same
# post tag commit, so its own recorded hash would be stale by construction.
# Drift in them is caught by continuous integration on every push instead.
SNAP = os.path.join(ROOT, "docs", "archive_snapshot.tsv")
n += 1
if not os.path.exists(SNAP):
    fails.append("archive snapshot: docs/archive_snapshot.tsv is absent, so the "
                 "deposit the manuscript cites cannot be checked against the "
                 "documents it is supposed to contain")
else:
    snap_tag = snap_doi = None
    recorded = []
    for raw in open(SNAP, encoding="utf-8"):
        line = raw.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts[0] == "tag":
            snap_tag = parts[1].strip()
        elif parts[0] == "version_doi":
            snap_doi = parts[1].strip()
        elif parts[0] == "file" and len(parts) >= 3:
            recorded.append((parts[1].strip(), parts[2].strip()))
    n += 1
    if snap_tag not in tags:
        fails.append(f"archive snapshot: the manuscript cites tag {tags!r} while "
                     f"the deposit on record is {snap_tag!r}")
    n += 1
    if not snap_doi or snap_doi not in text:
        fails.append(f"archive snapshot: the manuscript does not cite {snap_doi!r}, "
                     f"which is the DOI of the deposit on record")
    n += 1
    if not recorded:
        fails.append("archive snapshot: the record names no document, so it "
                     "asserts nothing")
    for rel, want in recorded:
        n += 1
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            fails.append(f"archive snapshot: {rel} is recorded as part of the "
                         f"deposit and is missing from the tree")
            continue
        body = open(path, encoding="utf-8").read()
        body = re.sub(r"10\.5281/zenodo\.\d+", "<DOI>", body)
        body = re.sub(r"\bv\d+\.\d+(\.\d+)?\b", "<TAG>", body)
        if hashlib.sha256(body.encode("utf-8")).hexdigest() != want:
            fails.append(f"archive snapshot: {rel} has changed since the deposit "
                         f"the manuscript cites was cut, so the cited DOI points "
                         f"at a different document. Cut a release and run "
                         f"make archive-snapshot TAG=... DOI=...")

# -- journal limits ---------------------------------------------------------
# Basis: Introduction through the start of the Declarations block, excluding
# table rows and display item captions. The limit excludes the abstract, the key
# messages, the declarations, the references, tables, figures and supplementary
# material, and this region is exactly what is left. No offset is applied: an
# earlier basis ran to the data availability statement and subtracted a
# calibration constant, which counted part of the declarations as main text and
# made the reported figure depend on the order of the declarations. Counting the
# region itself cannot understate, which is the direction that matters.
body = re.sub(r"^\|.*$", "", main_text, flags=re.M)
body = re.sub(r"^\*\*(Table|Figure) \d.*$", "", body, flags=re.M)
mw = len([w for w in re.split(r"\s+", re.sub(r"[*_#]", "", body)) if w.strip()])
n += 1
if LIMITS["main_text_words"] and mw > LIMITS["main_text_words"]:
    fails.append(f"journal limit: main text {mw} words exceeds the "
                 f"{LIMITS['main_text_words']} word limit")
aw = abstract_region
for h in ("**Abstract**", "**Background.**", "**Methods.**", "**Results.**", "**Conclusions.**"):
    aw = aw.replace(h, "")
aw = len([w for w in aw.split() if w.strip()])
n += 1
if LIMITS["abstract_words"] and aw > LIMITS["abstract_words"]:
    fails.append(f"journal limit: abstract {aw} words exceeds the "
                 f"{LIMITS['abstract_words']} word limit")
n += 1
title = text.splitlines()[0].strip().strip("*")
if LIMITS["title_characters"] and len(title) > LIMITS["title_characters"]:
    fails.append(f"journal limit: the title is {len(title)} characters, over the "
                 f"{LIMITS['title_characters']} character limit")

# -- share claims must carry a bound ---------------------------------------
# The first editorial comment was that the conclusions overstated the evidence,
# and the title was corrected on 19 October. The same claim survived in two
# Discussion sentences, "much of the genetic evidence" and "much of its apparent
# cancer risk", each asserting a majority share that the committed range of 32.5
# to 50.2 per cent does not support at its lower end. A denylist of phrases would
# have caught one wording and missed the other, so this is a property: a sentence
# that asserts a magnitude about the channel must carry a bound in the same
# sentence. Widen MAGNITUDE or SUBJECT, never add an exception for a wording.
MAGNITUDE = ("much", "most", "majority", "largely", "predominantly", "mostly",
             "the bulk")
SUBJECT = ("psa channel", "genetic signal", "genetic evidence", "that channel",
           "apparent cancer risk", "apparent risk", "apparent effect")
BOUNDING = ("part of", "a part", "between", "per cent", "%", "upper bound",
            "how much", "range")
for _sent in re.split(r"(?<=[.!?])\s+", main_text):
    _low = _sent.lower()
    if not any(re.search(r"\b%s\b" % _m, _low) for _m in MAGNITUDE):
        continue
    if not any(_s in _low for _s in SUBJECT):
        continue
    n += 1
    if not any(_b in _low for _b in BOUNDING):
        fails.append("share claim: a magnitude is asserted about the PSA channel "
                     "with no bound in the same sentence: %r"
                     % " ".join(_sent.split())[:140])

# -- figure conformance -----------------------------------------------------
# The Endocrine Society Digital Art Guidelines fix the publication size: one
# column 3.25 in, two column 5.8 in, height at most 9.25 in, submitted at final
# size and not scaled. Every figure in this repository was over the two column
# width until 19 October, the widest by 76 per cent, which would have had
# production scale an 8 point label to 4.8 point. figures/_style.py refuses to
# write an oversized figure, but that only runs when the figures are rebuilt.
# This reads what is actually committed.
#
# The PNG header is parsed directly rather than through an imaging library, so
# the check adds no dependency: IHDR carries the pixel dimensions and pHYs the
# pixels per unit, with unit 1 meaning metres.
FIG_MAX_W_IN = 5.8
FIG_MAX_H_IN = 9.25
FIG_MIN_DPI = 600


def _png_geometry(path):
    with open(path, "rb") as fh:
        if fh.read(8) != b"\x89PNG\r\n\x1a\n":
            return None
        px = py = wid = hgt = None
        while True:
            head = fh.read(8)
            if len(head) < 8:
                break
            ln = int.from_bytes(head[:4], "big")
            typ = head[4:8]
            data = fh.read(ln)
            fh.read(4)
            if typ == b"IHDR":
                wid = int.from_bytes(data[0:4], "big")
                hgt = int.from_bytes(data[4:8], "big")
            elif typ == b"pHYs" and data[8] == 1:
                px = int.from_bytes(data[0:4], "big") * 0.0254
                py = int.from_bytes(data[4:8], "big") * 0.0254
            elif typ == b"IDAT":
                break
        if None in (wid, hgt, px, py) or px <= 0 or py <= 0:
            return None
        return wid / px, hgt / py, px


_figdir = os.path.join(ROOT, "figures")
if os.path.isdir(_figdir):
    for _f in sorted(os.listdir(_figdir)):
        if not _f.endswith(".png"):
            continue
        n += 1
        _geom = _png_geometry(os.path.join(_figdir, _f))
        if _geom is None:
            fails.append("figure: %s carries no resolution, so its publication "
                         "size cannot be checked" % _f)
            continue
        _w, _h, _dpi = _geom
        if _w > FIG_MAX_W_IN + 0.01:
            fails.append("figure: %s is %.2f in wide at its stated resolution, "
                         "over the %.2f in two column maximum" % (_f, _w, FIG_MAX_W_IN))
        elif _h > FIG_MAX_H_IN + 0.01:
            fails.append("figure: %s is %.2f in tall, over the %.2f in maximum"
                         % (_f, _h, FIG_MAX_H_IN))
        elif _dpi < FIG_MIN_DPI - 1:
            fails.append("figure: %s is %d dpi, under the %d dpi the guidelines "
                         "give for a figure carrying line art and text"
                         % (_f, round(_dpi), FIG_MIN_DPI))

# -- title page concordance -------------------------------------------------
# The title page states six counts and closes by claiming they are asserted
# mechanically and cannot drift. Nothing read that page until this gate existed,
# and three of the six had drifted: the title, the main text length and the
# abstract length were those of a superseded draft. The counts are recomputed
# here from the manuscript, so the claim the page makes about itself is the
# claim this module enforces.
TITLE_PAGE = os.path.join(ROOT, "manuscript", "title_page.md")
if os.path.isfile(TITLE_PAGE):
    tp = open(TITLE_PAGE, encoding="utf-8").read()
    n += 1
    _m = re.search(r"^\*\*Title\.\*\*\s*(.+?)\s*$", tp, re.M)
    if not _m:
        fails.append("title page: states no title")
    elif _m.group(1) != title:
        fails.append("title page: the stated title is not the manuscript title")
    _tables = len(set(re.findall(r"\*\*Table (\d+)[.,]", text)))
    _figures = len(set(re.findall(r"\*\*Figure (\d+)[.,]", text)))
    _refblock = re.search(r"\*\*References\*\*(.*)$", text, re.S)
    _refs = (len(re.findall(r"^\s*\d+\.\s", _refblock.group(1), re.M))
             if _refblock else 0)
    for _label, _pat, _computed in (
            ("title characters", r"Title ([\d,]+) characters", len(title)),
            ("main text words", r"Main text ([\d,]+) words", mw),
            ("abstract words", r"Abstract ([\d,]+) words", aw),
            ("tables", r"([\d,]+) tables", _tables),
            ("figures", r"([\d,]+) figures", _figures),
            ("references", r"([\d,]+) references", _refs)):
        n += 1
        _c = re.search(_pat, tp)
        if not _c:
            fails.append("title page: states no %s count" % _label)
        elif int(_c.group(1).replace(",", "")) != _computed:
            fails.append("title page: states %s for %s, the manuscript gives %d"
                         % (_c.group(1), _label, _computed))

# -- submission gate: unresolved placeholders -------------------------------
PLACEHOLDERS = re.findall(r"\[[A-Z][A-Z0-9_]{6,}\]", text)
if SUBMISSION:
    n += 1
    if PLACEHOLDERS:
        fails.append("submission gate: unresolved placeholders remain: "
                     + ", ".join(sorted(set(PLACEHOLDERS))))

print("main text %d words%s, abstract %d of %s words, title %d of %s characters."
      % (mw,
         "" if not LIMITS["main_text_words"] else " of %d" % LIMITS["main_text_words"],
         aw, LIMITS["abstract_words"] or "no", len(title),
         LIMITS["title_characters"] or "no"))
if PLACEHOLDERS and not SUBMISSION:
    print("PENDING before submission: " + ", ".join(sorted(set(PLACEHOLDERS))))

print(f"{n} manuscript values checked against the committed matrix.")
if fails:
    print("\nDISCORDANCE:")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print("Abstract and Tables 2, 3 and 4 are concordant with the committed matrix.")
print("Claim consistency, count qualification, supplementary citation and "
      "journal limits all pass.")
