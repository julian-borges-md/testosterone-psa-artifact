#!/usr/bin/env python3
"""Between cohort comparison of the testosterone effect on prostate cancer.

The two cancer sources differ in how cases are found. PRACTICAL is a case
control consortium whose cases are largely PSA screen detected, with a case
fraction of 0.564. FinnGen is registry ascertained in a population without an
organised PSA screening programme, with a case fraction of 0.115.

If the apparent effect of testosterone on diagnosed cancer runs through PSA
driven detection, it must be larger where detection is PSA driven. That is a
cohort level prediction of the same mechanism the mediation analysis tests,
and it is independent of it. This module computes the comparison.

It is a post hoc consistency check and is reported as one. It was not among
the prespecified analyses, and nothing in the decision rules depends on it.

Writes results/stage7_matrix/between_cohort.tsv and prints the summary.
"""
import csv
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MATRIX = os.path.join(ROOT, "results", "stage7_matrix", "estimator_matrix.tsv")
OUT = os.path.join(ROOT, "results", "stage7_matrix", "between_cohort.tsv")

SCREENED, REGISTRY = "biot->prac", "biot->fg"


def norm_sf(z):
    return 0.5 * math.erfc(abs(z) / math.sqrt(2))


def load():
    by = {}
    with open(MATRIX, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by.setdefault(row["analysis"], {})[row["estimator"]] = (
                float(row["beta"]), float(row["se"]))
    return by


def main():
    by = load()
    for a in (SCREENED, REGISTRY):
        if a not in by:
            raise SystemExit("analysis %s absent from the committed matrix" % a)

    rows = []
    for est in by[SCREENED]:
        if est not in by[REGISTRY]:
            continue
        b1, s1 = by[SCREENED][est]
        b2, s2 = by[REGISTRY][est]
        d = b1 - b2
        sd = math.sqrt(s1 ** 2 + s2 ** 2)
        z = d / sd
        p = 2 * norm_sf(z)
        w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
        pooled = (w1 * b1 + w2 * b2) / (w1 + w2)
        q = w1 * (b1 - pooled) ** 2 + w2 * (b2 - pooled) ** 2
        i2 = max(0.0, (q - 1) / q) * 100 if q > 0 else 0.0
        rows.append({
            "estimator": est,
            "or_screen_detected": round(math.exp(b1), 4),
            "or_registry": round(math.exp(b2), 4),
            "difference_log_odds": round(d, 5),
            "se_difference": round(sd, 5),
            "z": round(z, 3),
            "p": round(p, 4),
            "cochran_q": round(q, 3),
            "i2_percent": round(i2, 1),
            "direction_as_predicted": "yes" if d > 0 else "no",
        })

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    predicted = sum(1 for r in rows if r["direction_as_predicted"] == "yes")
    sig05 = sum(1 for r in rows if r["p"] < 0.05)
    sig10 = sum(1 for r in rows if r["p"] < 0.10)
    print("between cohort comparison, bioavailable testosterone on prostate cancer")
    print("  screen detected source : %s" % SCREENED)
    print("  registry source        : %s" % REGISTRY)
    print()
    print("  %-28s %9s %9s %8s %7s %7s" %
          ("estimator", "screened", "registry", "diff", "p", "I2"))
    for r in rows:
        print("  %-28s %9.3f %9.3f %+8.4f %7.3f %6.1f%%" %
              (r["estimator"], r["or_screen_detected"], r["or_registry"],
               r["difference_log_odds"], r["p"], r["i2_percent"]))
    print()
    print("  larger in the screen detected source : %d of %d estimators" % (predicted, n))
    print("  difference significant at p<0.05     : %d of %d" % (sig05, n))
    print("  difference significant at p<0.10     : %d of %d" % (sig10, n))
    print("  range of p                           : %.3f to %.3f"
          % (min(r["p"] for r in rows), max(r["p"] for r in rows)))
    print("  range of between cohort I squared    : %.1f%% to %.1f%%"
          % (min(r["i2_percent"] for r in rows), max(r["i2_percent"] for r in rows)))
    print()
    print("  wrote %s" % os.path.relpath(OUT, ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
