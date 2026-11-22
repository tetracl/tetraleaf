#!/usr/bin/env python3
"""Runs a command and appends its peak resident set to a file (Helmholtz CH4.4a), in place
of GNU time, which Kathleen's compute nodes lack:

    python3 tools/fmm3d/maxrss.py FILE LABEL -- COMMAND [ARGUMENT...]

After COMMAND exits, it appends the line `maxrss: LABEL kb=<n>` to FILE, n the largest
resident set in KiB of the command and the descendants it waited for (getrusage's
RUSAGE_CHILDREN; KiB on Linux, bytes on macOS, converted), and exits with the command's
status. `tools/fmm3d/run.sh helmholtz` runs FMM3D's driver under it with LABEL
`token=<eps>`; `tetraleaf_validate::helmholtz_fmm3d::maxrss` reads the line. Python 3.9
standard library only.
"""

import resource
import subprocess
import sys


def main():
    if len(sys.argv) < 5 or sys.argv[3] != "--":
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    path, label, command = sys.argv[1], sys.argv[2], sys.argv[4:]
    status = subprocess.call(command)
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if sys.platform == "darwin":
        peak //= 1024
    with open(path, "a", encoding="utf-8") as out:
        out.write(f"maxrss: {label} kb={peak}\n")
    sys.exit(status if status >= 0 else 128 - status)


if __name__ == "__main__":
    main()
