#!/usr/bin/env python3
"""Lint the Fortran lines changed relative to a base revision with Fortitude.

Only lines added or modified since the merge base of HEAD and --base are
reported, so existing code is not reformatted as a side effect of unrelated
changes.
New files, including untracked ones, are therefore linted in full.
Rules, severities and excluded paths come from fortitude.toml in the
repository root.
Error-severity diagnostics make the script exit with a non-zero status;
warnings are reported but never change the exit status.

Run from anywhere in the repository, for example:
  python3 .github/scripts/fortitude_lint.py
  python3 .github/scripts/fortitude_lint.py --base upstream/main
  python3 .github/scripts/fortitude_lint.py --whole-files
Arguments after "--" are passed to "fortitude check", for example:
  python3 .github/scripts/fortitude_lint.py -- --output-format full

Fortitude 0.10 or newer is required:
  pip install 'fortitude-lint==0.10.*'
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys

FORTRAN_SUFFIXES = (".F", ".f", ".F90", ".f90", ".pf")
# "@@ -a[,b] +c[,d] @@": the new-file range starts at line c and spans d lines (default 1).
HUNK_RE = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")


def git(*args: str) -> str:
    """Run a git command and return its stdout, exiting on failure."""
    result = subprocess.run(["git", *args], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        sys.exit("git {} failed:\n{}".format(" ".join(args), result.stderr))
    return result.stdout


def default_base() -> str:
    """Prefer the remote-tracking main branch, which is current in CI and fresh clones."""
    for candidate in ("origin/main", "main"):
        check = subprocess.run(["git", "rev-parse", "--verify", "--quiet", candidate + "^{commit}"],
                               capture_output=True, check=False)
        if check.returncode == 0:
            return candidate
    sys.exit("Could not find origin/main or main; pass --base explicitly.")


def is_fortran(path: str) -> bool:
    return path.endswith(FORTRAN_SUFFIXES)


def changed_lines(merge_base: str) -> dict[str, list[list[int]]]:
    """Map each changed Fortran file to the inclusive line ranges added in the working tree."""
    # Compare against the working tree so uncommitted edits are linted too.
    # --diff-filter=d drops deleted files; -M keeps renamed files to their edited lines only.
    diff = git("diff", "-U0", "--no-color", "--no-ext-diff", "-M", "--diff-filter=d", merge_base, "--")
    ranges: dict[str, list[list[int]]] = {}
    path = None
    for line in diff.splitlines():
        if line.startswith("+++ "):
            target = line[4:]
            path = target[2:] if target.startswith("b/") else None
            if path is not None and not is_fortran(path):
                path = None
        elif path is not None:
            match = HUNK_RE.match(line)
            if match:
                start = int(match.group(1))
                count = int(match.group(2)) if match.group(2) is not None else 1
                # A pure deletion has count 0 and adds no lines to lint.
                if count > 0:
                    ranges.setdefault(path, []).append([start, start + count - 1])
    return ranges


def untracked_fortran() -> list[str]:
    """New files not yet added to git; they are linted in full."""
    return [path for path in git("ls-files", "--others", "--exclude-standard").splitlines()
            if is_fortran(path)]


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--base", default=None,
                        help="revision to compare against (default: origin/main, else main)")
    parser.add_argument("--whole-files", action="store_true",
                        help="lint every line of the changed files, not just the changed lines")
    parser.add_argument("--fortitude", default="fortitude",
                        help="fortitude executable (default: fortitude on PATH)")
    parser.add_argument("fortitude_args", nargs="*",
                        help="extra arguments for 'fortitude check', given after '--'")
    args = parser.parse_args(argv)

    fortitude = shutil.which(args.fortitude)
    if fortitude is None:
        sys.exit("Cannot find '{}'; install it with: pip install 'fortitude-lint==0.10.*'"
                 .format(args.fortitude))

    root = git("rev-parse", "--show-toplevel").strip()
    base = args.base or default_base()
    merge_base = git("merge-base", base, "HEAD").strip()

    # Paths from git are relative to the repository root, so run everything from there.
    os.chdir(root)

    ranges = changed_lines(merge_base)
    files = sorted(set(ranges) | set(untracked_fortran()))
    if not files:
        print("No Fortran files changed relative to {} ({}).".format(base, merge_base[:12]))
        return 0

    print("Linting {} Fortran file(s) changed relative to {} ({}).".format(
        len(files), base, merge_base[:12]), flush=True)

    command = [fortitude, "check", "--config", "fortitude.toml"]
    if not args.whole_files:
        # Files listed without "lines" (new or untracked) are linted in full.
        line_filter = [{"name": path, "lines": ranges[path]} if path in ranges else {"name": path}
                       for path in files]
        command += ["--line-filter", json.dumps(line_filter, separators=(",", ":"))]
    command += args.fortitude_args
    command += files
    result = subprocess.run(command, stdout=subprocess.PIPE, text=True, check=False)
    # Fortitude reports absolute paths, but GitHub annotations need paths relative to the repository root.
    output = result.stdout
    for prefix in {root, os.path.realpath(root)}:
        output = output.replace(prefix + os.sep, "")
    sys.stdout.write(output)
    return result.returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
