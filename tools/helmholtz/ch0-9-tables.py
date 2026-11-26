#!/usr/bin/env python3
"""tools/helmholtz/ch0-9-tables.py: the time tables of Helmholtz CH0.9 (spike S9,
docs/helmholtz/tasks/CH0.9.md) from one machine's `times-<variant>-<grid>.md` files
(spikes/fft-plans/run.sh writes them into spikes/fft-plans/results/<machine>/).

    tools/helmholtz/ch0-9-tables.py <results-dir>
    tools/helmholtz/ch0-9-tables.py selftest

Prints two tables, per degree p:
  times   the φ step per box, S → F plus F → R, in µs: the product's (its whole conversions
          times the φ share of their multiply-adds, CH4.1a's method) and each FFT back end on
          each grid (the scatter and gather included)
  ratios  FFTW and rustfft on the smooth grid over the product's φ step; FFTW's 2P rings over
          2P − 1, split over interleaved, the power-of-two and the product grid over the smooth
          one; rustfft's grids over its smooth one; and a whole conversion with the faster FFT
          on the smooth grid in place of the φ step, over the product's:
          ((1 − share) · product + FFT) / product

Python 3.9 standard library only; the tables are read with tools/helmholtz/common.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

FFT_COLUMNS = [
    ("fftw", "smooth", "FFTW smooth"),
    ("fftw-even", "smooth", "FFTW smooth 2P"),
    ("fftw-split", "smooth", "FFTW split smooth"),
    ("fftw", "pow2", "FFTW pow2"),
    ("fftw", "product", "FFTW product grid"),
    ("rustfft", "smooth", "rustfft smooth"),
    ("rustfft", "pow2", "rustfft pow2"),
    ("rustfft", "product", "rustfft product grid"),
]


def read(path):
    """{p: row} of one times file's table, or {} if the file is missing."""
    if not Path(path).exists():
        return {}
    found = common.tables(path)
    table = common.find_table(found, first="p")
    return {int(r["p"]): r for r in table.rows} if table else {}


def both(row):
    """S → F plus F → R µs of a row, or None."""
    if not row:
        return None
    a, b = common.number(row.get("S → F µs")), common.number(row.get("F → R µs"))
    return a + b if a is not None and b is not None else None


def ratio(a, b):
    return f"{a / b:.2f}" if a is not None and b else "–"


def us(x):
    return "–" if x is None else (f"{x:.0f}" if x >= 100 else f"{x:.1f}")


def report(directory):
    d = Path(directory)
    product = read(d / "times-product-product.md")
    ffts = {(v, g): read(d / f"times-{v}-{g}.md") for v, g, _ in FFT_COLUMNS}
    degrees = sorted(product)
    lines = []
    header = ["p", "N_φ product, smooth, pow2", "product whole", "φ share", "product φ"]
    header += [name for _, _, name in FFT_COLUMNS]
    lines.append("### φ step per box, S → F + F → R, µs\n")
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + " --- |" * len(header))
    for p in degrees:
        row = product[p]
        whole = both(row)
        share = common.number(row["φ share"])
        grids = [
            row["N_φ"],
            ffts[("fftw", "smooth")].get(p, {}).get("N_φ", "–"),
            ffts[("fftw", "pow2")].get(p, {}).get("N_φ", "–"),
        ]
        cells = [str(p), ", ".join(grids), us(whole), f"{share:.3f}", us(whole * share)]
        cells += [us(both(ffts[(v, g)].get(p))) for v, g, _ in FFT_COLUMNS]
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("\n### Ratios\n")
    header = [
        "p",
        "FFTW smooth ÷ product φ",
        "rustfft smooth ÷ product φ",
        "FFTW 2P ÷ 2P − 1",
        "FFTW split ÷ interleaved",
        "FFTW pow2 ÷ smooth",
        "FFTW product grid ÷ smooth",
        "rustfft pow2 ÷ smooth",
        "rustfft product grid ÷ smooth",
        "conversion with the faster FFT ÷ product",
    ]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + " --- |" * len(header))
    for p in degrees:
        whole = both(product[p])
        share = common.number(product[p]["φ share"])
        phi = whole * share
        t = {k: both(v.get(p)) for k, v in ffts.items()}
        fftw, rust = t[("fftw", "smooth")], t[("rustfft", "smooth")]
        best = min(x for x in (fftw, rust) if x is not None) if (fftw or rust) else None
        cells = [
            str(p),
            ratio(fftw, phi),
            ratio(rust, phi),
            ratio(t[("fftw-even", "smooth")], fftw),
            ratio(t[("fftw-split", "smooth")], fftw),
            ratio(t[("fftw", "pow2")], fftw),
            ratio(t[("fftw", "product")], fftw),
            ratio(t[("rustfft", "pow2")], rust),
            ratio(t[("rustfft", "product")], rust),
            ratio(None if best is None else (1 - share) * whole + best, whole),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def selftest():
    import tempfile

    with tempfile.TemporaryDirectory() as t:
        d = Path(t)
        (d / "times-product-product.md").write_text(
            "## times product product\n\nMachine: x\n\n"
            "| p | N_φ | boxes | S → F µs | F → R µs | φ share | φ S → F µs | φ F → R µs |\n"
            "| --- | --- | --- | --- | --- | --- | --- | --- |\n"
            "| 33 | 136 | 16 | 90.0 | 110.0 | 0.800 | 72.0 | 88.0 |\n",
            encoding="utf-8",
        )
        (d / "times-fftw-smooth.md").write_text(
            "## times fftw smooth\n\n| p | N_φ | rings | boxes | S → F µs | F → R µs |\n"
            "| --- | --- | --- | --- | --- | --- |\n| 33 | 140 | 67 | 16 | 20.0 | 20.0 |\n",
            encoding="utf-8",
        )
        (d / "times-rustfft-smooth.md").write_text(
            "## times rustfft smooth\n\n| p | N_φ | rings | boxes | S → F µs | F → R µs |\n"
            "| --- | --- | --- | --- | --- | --- |\n| 33 | 140 | 67 | 16 | 30.0 | 34.0 |\n",
            encoding="utf-8",
        )
        out = report(d)
    rows = [common.cells(l) for l in out.splitlines() if l.startswith("| 33 ")]
    assert rows[0][:5] == ["33", "136, 140, –", "200", "0.800", "160"], rows[0]
    assert rows[0][5] == "40.0" and rows[0][10] == "64.0", rows[0]
    # FFTW 40 / 160; rustfft 64 / 160; conversion (0.2 · 200 + 40) / 200.
    assert rows[1][1:3] == ["0.25", "0.40"] and rows[1][-1] == "0.40", rows[1]
    assert rows[1][3] == "–", rows[1]
    print("selftest passed")


if __name__ == "__main__":
    if sys.argv[1:] == ["selftest"]:
        selftest()
    elif len(sys.argv) == 2:
        print(report(sys.argv[1]))
    else:
        sys.exit(__doc__)
