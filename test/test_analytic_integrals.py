"""Compare each PDF's analytic integral with RooFit's numerical integral of evaluate()."""

import math
import os
import random
import sys

import ROOT

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")


def _report_mismatches(pdf, x, ranges, rtol):
    """
    Yields
      '[range] analytic=.. numeric=..'
    wherever the two integrals over x disagree.
    """
    obs = ROOT.RooArgSet(x)
    for rng in ranges:
        analytic, numeric = (_calculateIntegral(pdf, obs, rng, force) for force in (False, True))
        if not (math.isfinite(analytic) and math.isclose(analytic, numeric, rel_tol=rtol)):
            yield f"[{rng or 'full'}] analytic={analytic:.6g} numeric={numeric:.6g}"


def _calculateIntegral(pdf, obs, rng, force_numeric):
    pdf.forceNumInt(force_numeric)
    val = (pdf.createIntegral(obs, rng) if rng else pdf.createIntegral(obs)).getVal()
    pdf.forceNumInt(False)
    return val


def _has_analytic(pdf, x):
    return pdf.getAnalyticalIntegral(ROOT.RooArgSet(x), ROOT.RooArgSet()) != 0


def check_integral(cls, x, params, ranges=(None,), rtol=1e-4):
    """
    Build ROOT.<cls>(name, title, x, *params) and return where the analytic disagrees with the numeric integral.

    `ranges` holds named ranges already defined on `x` to check. `None` corresponds to the full range.
    Any class without an analytic integral over x is skipped.
    """
    pdf = getattr(ROOT, cls)(cls, cls, x, *params)
    return [f"{cls} {m}" for m in _report_mismatches(pdf, x, ranges, rtol)] if _has_analytic(pdf, x) else []


def check_integrals(cases, **kwargs):
    """
    Check the integrals for each case in `cases`.

    `cases` is a list of (cls, x, params) tuples. Returns all failures.
    """
    return [f for case in cases for f in check_integral(*case, **kwargs)]


def discover_pdfs(lib="libHiggsAnalysisCombinedLimit"):
    """
    Discover all RooAbsPdf classes in `lib` with a (name, title, RooAbsReal&...) constructor.

    Return {cls: number of RooAbsReal& args} for every RooAbsPdf in `lib`.
    Assumes the first RooAbsReal& is the observable."""
    with open(os.path.splitext(ROOT.gSystem.DynamicPathName(lib))[0] + ".rootmap") as f:
        classes = sorted({line.split()[1] for line in f if line.startswith("class ")})
    found = {}
    for cls in classes:
        tcl = ROOT.TClass.GetClass(cls)
        if not tcl or not tcl.InheritsFrom("RooAbsPdf") or tcl.Property() & ROOT.kIsAbstract:
            continue
        for meth in tcl.GetListOfMethods():
            types = [a.GetFullTypeName() for a in meth.GetListOfMethodArgs()]
            if meth.GetName() == cls and len(types) > 2 and types[:2] == ["const char*"] * 2 and set(types[2:]) == {"RooAbsReal&"}:
                found[cls] = len(types) - 2
                break
    return found


def _in_domain(pdf, x, n=50):
    """True if evaluate() is finite and non-negative across x's range, and not identically zero."""
    lo, hi = x.getMin(), x.getMax()
    vals = []
    for i in range(n):
        x.setVal(lo + (i + 0.5) * (hi - lo) / n)
        vals.append(pdf.getVal())
    return all(math.isfinite(v) and v >= 0 for v in vals) and any(vals)


def auto_check(n_points=20, seed=1, rtol=1e-4, param_range=(-5, 5)):
    """Check every discovered PDF at random parameter points inside its domain.

    Returns {cls: status}, where status is a mismatch, "ok (k points)", or a reason it was skipped.
    """
    rng = random.Random(seed)
    x = ROOT.RooRealVar("x", "x", 5, 1, 10)
    x.setRange("sub", 3, 6)
    results = {}
    for cls, n in discover_pdfs().items():
        params = [ROOT.RooRealVar(f"{cls}_p{i}", "", 1, -1e6, 1e6) for i in range(n - 1)]
        try:
            pdf = getattr(ROOT, cls)(cls, cls, x, *params)
        except Exception as e:  # noqa: BLE001 - any constructor failure means skip
            results[cls] = f"skipped: constructor failed ({str(e).splitlines()[-1].strip()[:60]})"
            continue
        if not _has_analytic(pdf, x):
            results[cls] = "skipped: no analytic integral"
            continue
        valid = 0
        for _ in range(n_points):
            for p in params:
                p.setVal(rng.uniform(*param_range))
            if not _in_domain(pdf, x):
                continue
            valid += 1
            if mismatch := next(_report_mismatches(pdf, x, (None, "sub"), rtol), None):
                results[cls] = f"MISMATCH {mismatch} at params=({', '.join(f'{p.getVal():.3g}' for p in params)})"
                break
        else:
            results[cls] = f"ok ({valid} points)" if valid else "skipped: no in-domain point found"
    return results


def regression_cases():
    """
    Return (cls, x, params) cases at the parameter points where analytic integrals used to be wrong.

    Each x has a "sub" range defined.
    """
    # (cls, x range, sub range, params)
    points = [
        ("Triangle", (0, 10), (4, 6), (-2, 5, 3)),
        ("Triangle", (0, 10), (4, 6), (-4, 2, 3)),  # cut off at x=0
        ("RooExpPdf", (110, 200), (120, 160), (-0.05, 91.2)),
        ("RooExpPdf", (110, 200), (120, 160), (0, 91.2)),
        ("RooSumTwoPowerLawPdf", (110, 200), (120, 160), (-2.5, -4, 0.3, 91.2)),
        ("RooSumTwoPowerLawPdf", (110, 200), (120, 160), (-1, -2, 0.5, 91.2)),  # log branches
        ("RooSumTwoPowerLawPdf", (110, 200), (120, 160), (-3, -1, 0.8, 91.2)),
        ("RooAnaExpNPdf", (500, 3000), (800, 1500), (-0.003, 0.5)),
        ("RooAnaExpNPdf", (500, 3000), (800, 1500), (-0.003, 100)),  # used to give NaN
        ("RooAnaExpNPdf", (500, 3000), (800, 1500), (-0.001, -50)),
        ("RooErfExpPdf", (1, 10), (3, 6), (-0.5, 3, 2)),
        ("RooErfExpPdf", (1, 10), (3, 6), (-4.75, 0.414, 4.39)),  # large |c| * width
        ("RooErfExpPdf", (1, 10), (3, 6), (4.6, 1.45, 3.27)),
    ]
    cases = []
    for i, (cls, (lo, hi), sub, values) in enumerate(points):
        x = ROOT.RooRealVar(f"x{i}", "x", lo, lo, hi)
        x.setRange("sub", *sub)
        params = [ROOT.RooRealVar(f"{cls}{i}_p{j}", "", v, -1e6, 1e6) for j, v in enumerate(values)]
        cases.append((cls, x, params))
    return cases


if __name__ == "__main__":
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    ROOT.RooAbsReal.setEvalErrorLoggingMode(ROOT.RooAbsReal.Ignore)
    # The default precision isn't tight enough for a reference at rtol=1e-6
    ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-12)
    ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-10)
    # Negative parameters put classes such as the Crystal Balls outside the domain
    # of their analytic formulas, so sample moderate positive values only.
    results = auto_check(param_range=(0.5, 3))
    for cls, status in sorted(results.items(), key=lambda kv: (not kv[1].startswith("MISMATCH"), kv[1], kv[0])):
        print(f"{cls:32s} {status}")
    failures = sorted(cls for cls, status in results.items() if status.startswith("MISMATCH"))
    regressions = check_integrals(regression_cases(), ranges=(None, "sub"), rtol=1e-6)
    for r in regressions:
        print(f"REGRESSION {r}")
    if failures or regressions:
        sys.exit(f"Analytic and numeric integrals disagree for: {', '.join(sorted(set(failures) | {r.split()[0] for r in regressions}))}")
    print("All checks passed")
