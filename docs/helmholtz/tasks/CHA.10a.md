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
- `tools/helmholtz/common.py` (new; as built, after the self-review): `cells(line)`;
  `number(text, default=None)`; `eps_token(text)` (the `1e-6` of `ε = 1e-6, …`); `bold(text,
  over)`; `gmean(values)` (positive finite values, the rest left out, NaN for none); `Table` and
  `tables(path)` (every Markdown table of an output in order, under its `## ` heading, rows as
  dicts, a repeated header name numbered `W (2)`); `find_table(found, heading, first)`;
  `accuracy_rows(path)` (`helmholtz_accuracy`'s run rows, the key "ε" added); `build_rows(path)`
  (the "Build" table's rows in order, "ε" added: one per ε in `helmholtz_fmm3d compare`, one per
  run in `helmholtz_accuracy`, where runs share an ε); `BoxesRun`, `boxes_runs(path)` and
  `boxes_run(path)` (`helmholtz_boxes`' run rows); `CompareRow` and `compare_rows(path)` (every
  row of `helmholtz_fmm3d compare`'s first table, ours and FMM3D's variants, with the cells as
  printed); `selftest()`, run as `python3 tools/helmholtz/common.py selftest`. Python 3.9
  standard library.
- `tools/helmholtz/ch4-8b-tables.py`: imports `common` (by its own directory, so `python3 -I`
  works) in place of its `table_cells`, `boxes_run`, `compare_rows`, `bold` and inline geometric
  means; a local `ours(path)` keeps our compare rows by ε; its `ratios`, `case_key`, tables, usage
  text and selftest stay (the selftest also runs `common.selftest()`).
- `.claude/skills/next-task/SKILL.md` step 5: the conditional "once CHA.10 has written it" goes;
  the plan template's "Changes" names the module for a report script.
- The old scripts are left as they are.

## Tests that prove it
- `python3 tools/helmholtz/common.py selftest`: each parser on lines of its format taken from the
  archived outputs (a `–` cell, a failed FMM3D call, a boxes row without errors, a missing table,
  a header with `|v|`, a repeated header name, two accuracy runs at one ε), `gmean` against a hand
  value and its left-out values, `number`'s fallback.
- Byte-identical output of `ch4-8b-tables.py` before (6bad3c2) and after, every mode over every
  archived directory it reads: Kathleen `ch4-8b` (`const`, `random` with and without CH4.4a's
  directory), `ch4-8b-cost` on both machines (`cost`), both machines' `h4-gate2` (`const`, W1–W3
  and W4), the `selftest` and the usage text: 34 outputs, `diff -r` empty, under Python 3.15 and
  3.9.6.
- `common.py` on the other two formats, against the archive (a check script outside the repo):
  `accuracy_rows` and `build_rows` on CH4.3c3's 52 `verify-matrix` outputs agree with
  `ch4-3c3-tables.py`'s `table_rows` (348 rows: ε, φ, ∇φ, checksum, degrees; the Build rows one per
  run in order); on all 686 `helmholtz_accuracy` outputs of the archive one Build row per run and
  no header name lost; `compare_rows` on CH4.4a's 26 outputs equals the raw cells (204 rows, ours
  and both FMM3D variants, the failed calls as None).

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
