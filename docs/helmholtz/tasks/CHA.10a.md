# CHA.10a The report scripts' shared module

Started 2026-11-25 on `tbetcke/h-cha-10a`. Roadmap: "CHA.10a The report scripts' shared module.
Added at the H4 gate's second run (the user on #262, 2026-11-25: "CHA.10a, then FFT, H2, H5
(Recommended)"; review 2026-11-25, K2, asked for at the H1 gate and the first H4 run).
`tools/helmholtz/common.py`: the cell formatting (`cells`, in fifteen `*-tables.py`), the parsers of
the harness examples' Markdown outputs (`helmholtz_boxes`, `helmholtz_accuracy`, `helmholtz_fmm3d
compare`) and the geometric means, with a `selftest`; the old scripts left as they are (their
reports' reproduction lines stay valid). *Done when:* the module exists with its selftest, one recent
script (`ch4-8b-tables.py`) uses it with byte-identical output, and `/next-task`'s planning says a
new report script imports it. Python standard library only. *Needs:* none." Proposal: none
(section HA: the roadmap's *Done when* with constraint 1 of §2).

## Goal
One standard-library module that new report scripts of `tools/helmholtz/` import for the Markdown
row splitting, the parsers of the three harness examples' outputs and the geometric mean, instead
of each defining its own.

## Acceptance (roadmap *Done when*, verbatim)
The module exists with its selftest, one recent script (`ch4-8b-tables.py`) uses it with
byte-identical output, and `/next-task`'s planning says a new report script imports it. Python
standard library only. With constraint 1 (the HA preamble): no product code changes.

## Read first
- Review `docs/reviews/2026-11-25-H4.md` K2 (checked on 2026-11-25 at 6bad3c2): `tools/helmholtz/`
  has 74 files; `def cells` (`[c.strip() for c in line.strip().strip("|").split("|")]`) in fourteen
  `*-tables.py` and as `table_cells` in `ch4-8b-tables.py` (`ch4-8a-tables.py`'s `cells` is a
  formatter of its own); a generic `tables(path)` ({heading: [row dicts]}) in fourteen scripts in
  six variants; `number`/`num` (float or a fallback) in eight, `gmean` in two, inline geometric
  means in several more.
- `tools/helmholtz/ch4-8b-tables.py`: `table_cells`, `boxes_run` (`helmholtz_boxes`' first run
  row), `compare_rows` (`helmholtz_fmm3d compare`'s rows of ours and its "Build" table), `bold`,
  the inline geometric means of `cost`, its `selftest`.
- The harness outputs as archived: `helmholtz_boxes` (`~/projects/fmm-archive/kathleen/runs/ch4-8b/
  const-*.md`), `helmholtz_fmm3d compare` (`…/ch4-8b-cost/cost-*.md`, `…/ch4-4a/helmholtz-*.md`),
  `helmholtz_accuracy` (`…/ch4-3c3/verify-matrix/*.md`: its first table, `| Workload | … | degrees
  | …`, with `ε = <eps>` in `degrees`, then `## <heading>` tables; its "Build" table keyed by
  `degrees`, compare's by `ε`).

## Changes
Scripts and documents only; no crate, manifest or product code.
- `tools/helmholtz/common.py` (new): `EPSILONS`; `cells(line)`; `number(text, default=None)`;
  `eps_token(text)` (the `1e-6` of `ε = 1e-6, …`); `tables(path)` (every Markdown table of an
  output by its `## ` heading, the first under `""`, rows as dicts); `accuracy_rows(path)`
  (`helmholtz_accuracy`'s main table, the ε added); `build_rows(path)` (the "Build" table by ε,
  both harnesses); `BoxesRun` and `boxes_runs(path)` / `boxes_run(path)` (`helmholtz_boxes`' run
  rows); `CompareRow` and `compare_rows(path, code="ours")` (`helmholtz_fmm3d compare`'s rows of
  one code by ε, and the "Build" table's build, evaluation and first evaluation by ε);
  `bold(text, over)`; `gmean(values)` (positive finite values; NaN for none); `selftest()` on
  lines of each format, run as `python3 tools/helmholtz/common.py selftest`.
- `tools/helmholtz/ch4-8b-tables.py`: imports `common` for the above in place of its
  `table_cells`, `boxes_run`, `compare_rows`, `bold` and inline geometric means; its own
  `ratios`, `case_key`, tables and selftest stay (the selftest also runs `common.selftest()`).
- `.claude/skills/next-task/SKILL.md` step 5: the conditional "once CHA.10 has written it" goes;
  the plan template's "Changes" names the module for a report script.
- The old scripts are left as they are.

## Tests that prove it
- `python3 tools/helmholtz/common.py selftest`: each parser on lines of its format (taken from the
  archived outputs, including a `–` cell, a missing table and a row whose cell holds `|`), `gmean`
  against a hand value, `number`'s fallback.
- Byte-identical output of `ch4-8b-tables.py` before and after, every mode over every archived
  directory it reads: Kathleen `ch4-8b` (`const`, `random` with and without CH4.4a's directory),
  `ch4-8b-cost` and locust's (`cost`), both machines' `h4-gate2` (`const`, W1–W3 and W4), the
  `selftest` and the usage text: 34 outputs captured at 6bad3c2, `diff -r` empty after.
- `common.py` on the other two formats: `accuracy_rows` and `build_rows` read every
  `helmholtz_accuracy` output of `…/ch4-3c3/verify-matrix/` and agree with
  `ch4-3c3-tables.py`'s `table_rows` (row count, ε and the φ, ∇φ cells), and `compare_rows`
  with `code="fmm3d"` reads CH4.5's outputs.

## Checks
`/verify` (scripts and documents: no cargo checks apply); `python3 -m py_compile` of both scripts;
the selftests; the byte-identity `diff -r`; `tools/roadmap/sync_issues.py --check`.

## Machine runs
None (the archived outputs are read locally).

## Constraints and risks
Constraint 1: no Laplace or Helmholtz code changes; the old scripts and their reports'
reproduction lines are untouched. The scripts run from any directory as
`tools/helmholtz/<name>.py` because Python puts the script's directory first on `sys.path`;
`common.py` is the only module there, `__pycache__/` is ignored by git.

## Pending runs
None.

## Open questions
None.
