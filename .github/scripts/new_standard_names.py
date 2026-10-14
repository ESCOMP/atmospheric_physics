#!/usr/bin/env python3
"""Report the CCPP standard names a pull request introduces or removes.

Standard names are collected from every tracked *.meta file
("standard_name = <name>" entries) and every tracked *.xml file
("<standard_name>" elements of the scheme namelist definitions).
Git submodules are not scanned.
A name is "new" when it occurs in the head revision but nowhere in the
merge base of head and base, which is exactly what the pull request diff
shows as added.
A name is "removed" when the opposite holds.
Adjacent removed/added standard_name lines in one diff hunk are reported
as a likely rename.
The units declared with each occurrence are reported too, so reviewers
can spot one name declared with different units.

Subcommands:
  report   compare two revisions and write report.json and report.md
           into --out-dir
  render   re-create the markdown from an existing report.json after
           validating it; used by the comment workflow, which must treat
           the downloaded report as untrusted input

Run locally from the repository root, for example:
  python3 .github/scripts/new_standard_names.py report \\
      --base origin/main --head HEAD --out-dir /tmp/stdnames
"""
from __future__ import annotations

import argparse
import difflib
import json
import os
import re
import subprocess
import sys
from collections import defaultdict

COMMENT_MARKER = "<!-- new-standard-names-report -->"
SCHEMA_VERSION = 1

# "standard_name = foo"; whole-line "#" comments never match because the
# key must come first, and a trailing comment is cut off by the class.
META_NAME_RE = re.compile(r"^\s*standard_name\s*=\s*([^\s#]+)")
XML_NAME_RE = re.compile(r"<standard_name>\s*([^<\s]+)\s*</standard_name>")
META_UNITS_RE = re.compile(r"^\s*units\s*=\s*([^#]*?)\s*(?:#.*)?$")
META_SECTION_RE = re.compile(r"^\s*\[")
XML_ENTRY_RE = re.compile(r"<entry\b.*?</entry>", re.DOTALL)
XML_UNITS_RE = re.compile(r"<units>\s*([^<]*?)\s*</units>")
SCANNED_SUFFIXES = (".meta", ".xml")

# Validation of report.json contents (see validate_report).
NAME_RE = re.compile(r"^[A-Za-z0-9_]{1,256}$")
PATH_RE = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_./-]{0,511}$")
REPO_RE = re.compile(r"^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
REF_RE = re.compile(r"^[A-Za-z0-9_./-]{1,256}$")
# Units end up inside markdown code spans, so backticks, pipes and angle
# brackets are never accepted.
UNITS_RE = re.compile(r"^[A-Za-z0-9 _.+*/^()%-]{1,64}$")
MAX_ENTRIES = 5000
MAX_OCCURRENCES = 10000
# A pull request comment may hold 65536 characters; keep the markdown under that.
MAX_MARKDOWN_CHARS = 60000


# ---------------------------------------------------------------------------
# git access
# ---------------------------------------------------------------------------
def git(*args: str, input_bytes: bytes | None = None) -> bytes:
    """Run a git command and return its stdout, raising on failure."""
    result = subprocess.run(
        ["git", *args], input=input_bytes, capture_output=True, check=False
    )
    if result.returncode != 0:
        sys.exit(
            "git {} failed:\n{}".format(
                " ".join(args), result.stderr.decode("utf-8", "replace")
            )
        )
    return result.stdout


def resolve_commit(rev: str) -> str:
    return git("rev-parse", "--verify", rev + "^{commit}").decode().strip()


def scanned_paths(rev: str) -> list[str]:
    """Tracked *.meta and *.xml blobs in <rev>; submodule gitlinks never match."""
    out = git("ls-tree", "-r", "-z", "--name-only", rev).decode("utf-8", "replace")
    return sorted(p for p in out.split("\0") if p.endswith(SCANNED_SUFFIXES))


def read_blobs(rev: str, paths: list[str]) -> dict[str, str]:
    """Fetch many blobs in one "git cat-file --batch" call.

    The batch output is "<sha> <type> <size>\\n<content>\\n" per object,
    or "<spec> missing\\n" for anything that does not resolve.
    """
    if not paths:
        return {}
    request = "".join("{}:{}\n".format(rev, p) for p in paths).encode()
    out = git("cat-file", "--batch", input_bytes=request)
    blobs: dict[str, str] = {}
    pos = 0
    for path in paths:
        newline = out.index(b"\n", pos)
        header = out[pos:newline].decode("utf-8", "replace").split()
        pos = newline + 1
        if header[-1] == "missing":
            continue
        obj_type, size = header[1], int(header[2])
        if obj_type == "blob":
            blobs[path] = out[pos:pos + size].decode("utf-8", "replace")
        pos += size + 1
    return blobs


# ---------------------------------------------------------------------------
# standard name extraction
# ---------------------------------------------------------------------------
def normalize_units(units: str | None) -> str | None:
    """Collapse internal whitespace so "m  s-1" and "m s-1" compare equal."""
    if units is None:
        return None
    units = " ".join(units.split())
    return units or None


def meta_names(text: str):
    """Yield (name, line_number, units) for each "[ var ]" entry of a .meta file.

    The units line may come before or after standard_name, so each entry is
    collected whole and emitted when the next section header begins.
    """
    pending: list[tuple[str, int]] = []
    units = None

    def flush():
        for name, number in pending:
            yield name, number, units

    for number, line in enumerate(text.splitlines(), 1):
        if META_SECTION_RE.match(line):
            yield from flush()
            pending, units = [], None
            continue
        match = META_NAME_RE.match(line)
        if match:
            pending.append((match.group(1), number))
            continue
        match = META_UNITS_RE.match(line)
        if match:
            units = normalize_units(match.group(1))
    yield from flush()


def xml_names(text: str):
    """Yield (name, line_number, units) for each <standard_name> in a namelist XML.

    Units come from the <units> element of the enclosing <entry>.
    """
    entries = [(m.start(), m.end(), m.group(0)) for m in XML_ENTRY_RE.finditer(text)]
    for match in XML_NAME_RE.finditer(text):
        units = None
        for start, end, body in entries:
            if start <= match.start() < end:
                found = XML_UNITS_RE.search(body)
                units = normalize_units(found.group(1)) if found else None
                break
        yield match.group(1), text.count("\n", 0, match.start()) + 1, units


def names_in_text(path: str, text: str):
    """Yield (name, line_number, units) for every standard name in one file."""
    if path.endswith(".meta"):
        yield from meta_names(text)
    else:
        yield from xml_names(text)


def scan_tree(rev: str) -> dict[str, list[tuple[str, int, str | None]]]:
    """Map every standard name in <rev> to its sorted (path, line, units) occurrences."""
    occurrences: dict[str, list[tuple[str, int, str | None]]] = defaultdict(list)
    for path, text in read_blobs(rev, scanned_paths(rev)).items():
        for name, line, units in names_in_text(path, text):
            occurrences[name].append((path, line, units))
    for locations in occurrences.values():
        locations.sort()
    return occurrences


def name_in_diff_line(content: str) -> str | None:
    match = META_NAME_RE.match(content) or XML_NAME_RE.search(content)
    return match.group(1) if match else None


def distinct_units(occurrences) -> list[str]:
    """Distinct declared units, most common first; "?" marks a missing units line."""
    counts: dict[str, int] = {}
    for occurrence in occurrences:
        units = occurrence[2] or "?"
        counts[units] = counts.get(units, 0) + 1
    return sorted(counts, key=lambda u: (-counts[u], u))


def detect_renames(merge_base: str, head: str, removed: set[str]) -> dict[str, str]:
    """Pair removed names with the names that take their place in a diff hunk.

    Returns {removed_name: replacement}; the replacement is either a new
    name (a rename) or one that already exists (a consolidation).
    Only hunks with the same number of removed and added standard_name
    lines are paired, in order, so "-standard_name = old" directly
    followed by "+standard_name = new" is caught while unrelated churn
    is left alone.
    """
    diff = git(
        "diff", "-U0", "--no-color", "--no-ext-diff",
        merge_base, head, "--", "*.meta", "*.xml",
    ).decode("utf-8", "replace")
    pairs: dict[str, str] = {}
    minus: list[str] = []
    plus: list[str] = []

    def flush_hunk():
        if minus and len(minus) == len(plus):
            for old, candidate in zip(minus, plus):
                if old in removed and candidate != old and old not in pairs:
                    pairs[old] = candidate
        minus.clear()
        plus.clear()

    for line in diff.splitlines():
        if line.startswith("@@") or line.startswith("diff --git"):
            flush_hunk()
        elif line.startswith("--- ") or line.startswith("+++ "):
            continue
        elif line.startswith("-"):
            name = name_in_diff_line(line[1:])
            if name:
                minus.append(name)
        elif line.startswith("+"):
            name = name_in_diff_line(line[1:])
            if name:
                plus.append(name)
    flush_hunk()
    return pairs


def similar_names(name: str, candidates: list[str]) -> list[str]:
    """Existing names that differ only slightly, to catch near-duplicates."""
    return difflib.get_close_matches(name, candidates, n=3, cutoff=0.9)


# ---------------------------------------------------------------------------
# report construction
# ---------------------------------------------------------------------------
def location(path: str, line: int, units: str | None = None) -> dict:
    loc = {"path": path, "line": line}
    if units is not None:
        loc["units"] = units
    return loc


def build_report(args: argparse.Namespace) -> dict:
    head = resolve_commit(args.head)
    base = resolve_commit(args.base)
    if args.merge_base:
        merge_base = resolve_commit(args.merge_base)
    else:
        merge_base = git("merge-base", base, head).decode().strip()

    head_names = scan_tree(head)
    merge_base_names = scan_tree(merge_base)
    base_tip_names = merge_base_names if base == merge_base else scan_tree(base)

    new_names = set(head_names) - set(merge_base_names)
    removed_names = set(merge_base_names) - set(head_names)
    replaced_by = detect_renames(merge_base, head, removed_names)
    renames = {new: old for old, new in replaced_by.items() if new in new_names}
    existing = sorted(merge_base_names)

    new_entries = []
    for name in sorted(new_names):
        occurrences = head_names[name]
        # Report the first .meta occurrence as where the name is introduced,
        # as that is where reviewers judge it; fall back to the XML.
        anchor = next((o for o in occurrences if o[0].endswith(".meta")), occurrences[0])
        old_name = renames.get(name)
        similar = similar_names(name, existing)
        new_entries.append({
            "name": name,
            "anchor": location(*anchor),
            "occurrences": [location(*o) for o in occurrences],
            "units": distinct_units(occurrences),
            "likely_rename_of": old_name,
            "likely_rename_of_units":
                distinct_units(merge_base_names[old_name]) if old_name else None,
            "on_base_tip": name in base_tip_names,
            "similar_existing": similar,
            "similar_existing_units":
                {other: distinct_units(merge_base_names[other]) for other in similar},
        })
    new_entries.sort(key=lambda e: (e["anchor"]["path"], e["anchor"]["line"]))

    removed_entries = [
        {
            "name": name,
            "occurrences": [location(*o) for o in merge_base_names[name]],
            "units": distinct_units(merge_base_names[name]),
            "likely_renamed_to": replaced_by.get(name),
        }
        for name in sorted(removed_names)
    ]

    base_ref = args.base_ref
    if not base_ref:
        base_ref = args.base[len("origin/"):] if args.base.startswith("origin/") else args.base
    return {
        "schema": SCHEMA_VERSION,
        "repository": args.repo,
        "pr": args.pr,
        "base_ref": base_ref,
        "base_sha": base,
        "merge_base_sha": merge_base,
        "head_sha": head,
        "new": new_entries,
        "removed": removed_entries,
    }


# ---------------------------------------------------------------------------
# validation of a report that may come from an untrusted workflow run
# ---------------------------------------------------------------------------
def check(condition: bool, what: str):
    if not condition:
        raise ValueError("invalid report: " + what)


def validate_location(loc) -> None:
    check(isinstance(loc, dict), "location is not an object")
    path, line = loc.get("path"), loc.get("line")
    check(isinstance(path, str) and PATH_RE.match(path) is not None, "bad path")
    check(".." not in path.split("/"), "path escapes the repository")
    check(isinstance(line, int) and 0 < line < 10_000_000, "bad line number")
    if "units" in loc:
        validate_units(loc["units"])


def validate_name(name) -> None:
    check(isinstance(name, str) and NAME_RE.match(name) is not None, "bad standard name")


def validate_units(units) -> None:
    check(isinstance(units, str) and UNITS_RE.match(units) is not None, "bad units")


def validate_units_list(units_list) -> None:
    """Units lists are optional so reports written before units existed still render."""
    if units_list is None:
        return
    check(isinstance(units_list, list) and 0 < len(units_list) <= 50, "bad units list")
    for units in units_list:
        validate_units(units)


def validate_report(report: dict) -> None:
    check(isinstance(report, dict), "not an object")
    check(report.get("schema") == SCHEMA_VERSION, "unsupported schema version")
    repo = report.get("repository")
    check(repo is None or (isinstance(repo, str) and REPO_RE.match(repo)), "bad repository")
    pr = report.get("pr")
    check(pr is None or (isinstance(pr, int) and 0 < pr < 10_000_000), "bad pr number")
    for key in ("base_sha", "merge_base_sha", "head_sha"):
        value = report.get(key)
        check(isinstance(value, str) and SHA_RE.match(value) is not None, "bad " + key)
    base_ref = report.get("base_ref")
    check(isinstance(base_ref, str) and REF_RE.match(base_ref) is not None, "bad base_ref")

    occurrences_seen = 0
    for key in ("new", "removed"):
        entries = report.get(key)
        check(isinstance(entries, list) and len(entries) <= MAX_ENTRIES, "bad " + key)
        for entry in entries:
            check(isinstance(entry, dict), "entry is not an object")
            validate_name(entry.get("name"))
            occurrences = entry.get("occurrences")
            check(isinstance(occurrences, list) and occurrences, "missing occurrences")
            occurrences_seen += len(occurrences)
            check(occurrences_seen <= MAX_OCCURRENCES, "too many occurrences")
            for loc in occurrences:
                validate_location(loc)
            validate_units_list(entry.get("units"))
            if key == "new":
                validate_location(entry.get("anchor"))
                check(entry.get("on_base_tip") in (True, False), "bad on_base_tip")
                rename = entry.get("likely_rename_of")
                if rename is not None:
                    validate_name(rename)
                validate_units_list(entry.get("likely_rename_of_units"))
                similar = entry.get("similar_existing")
                check(isinstance(similar, list) and len(similar) <= 3, "bad similar_existing")
                for other in similar:
                    validate_name(other)
                similar_units = entry.get("similar_existing_units")
                if similar_units is not None:
                    check(isinstance(similar_units, dict)
                          and set(similar_units) <= set(similar), "bad similar_existing_units")
                    for units_list in similar_units.values():
                        validate_units_list(units_list)
            else:
                rename = entry.get("likely_renamed_to")
                if rename is not None:
                    validate_name(rename)


# ---------------------------------------------------------------------------
# output: markdown (job summary and sticky comment)
# ---------------------------------------------------------------------------
def short(sha: str) -> str:
    return sha[:8]


def link(report: dict, loc: dict) -> str:
    text = "`{}:{}`".format(loc["path"], loc["line"])
    repo = report.get("repository")
    if repo:
        return "[{}](https://github.com/{}/blob/{}/{}#L{})".format(
            text, repo, report["head_sha"], loc["path"], loc["line"]
        )
    return text


def code_list(values: list[str]) -> str:
    return ", ".join("`{}`".format(v) for v in values)


def units_cell(entry: dict) -> str:
    """Units column text; blank for reports written before units were recorded."""
    return "<br>".join("`{}`".format(u) for u in entry.get("units") or [])


def units_differ(entry: dict) -> bool:
    return len(entry.get("units") or []) > 1


def location_list(report: dict, locations: list[dict], limit: int = 6,
                  with_units: bool = False) -> str:
    """Linked locations; with_units appends each one's units, used when they disagree."""
    def one(loc):
        text = link(report, loc)
        if with_units:
            text += " (`{}`)".format(loc.get("units", "?"))
        return text

    shown = [one(loc) for loc in locations[:limit]]
    if len(locations) > limit:
        shown.append("and {} more".format(len(locations) - limit))
    return "<br>".join(shown)


def render_markdown(report: dict, with_marker: bool = False) -> str:
    base_ref = report["base_ref"]
    new, removed = report["new"], report["removed"]
    lines = []
    if with_marker:
        lines.append(COMMENT_MARKER)
    lines.append("## Standard names introduced by this pull request")
    lines.append("")
    lines.append(
        "Non-blocking review aid, refreshed on every push. "
        "Names in every tracked `*.meta` and `*.xml` file (submodules excluded) "
        "at head {} are compared with `{}` at the merge base {}.".format(
            short(report["head_sha"]), base_ref, short(report["merge_base_sha"])
        )
    )
    lines.append("")
    if not new and not removed:
        lines.append("No standard names were added or removed.")
        return "\n".join(lines) + "\n"

    renamed = sum(1 for e in new if e["likely_rename_of"])
    summary = "**{} new standard name{}**".format(len(new), "" if len(new) == 1 else "s")
    if renamed:
        summary += " ({} likely rename{})".format(renamed, "" if renamed == 1 else "s")
    summary += " and **{} removed**.".format(len(removed))
    lines.append(summary)
    lines.append("")

    if new:
        lines.append("| New standard name | Units | Introduced in | Also seen in | Notes |")
        lines.append("|---|---|---|---|---|")
        table_start = len(lines)
        for entry in new:
            anchor = entry["anchor"]
            others = [o for o in entry["occurrences"] if o != anchor]
            differ = units_differ(entry)
            notes = []
            if differ:
                notes.append("**declared with different units**")
            if entry["likely_rename_of"]:
                note = "likely rename of `{}`".format(entry["likely_rename_of"])
                old_units = entry.get("likely_rename_of_units")
                if old_units and old_units != entry.get("units"):
                    note += ", which had units " + code_list(old_units)
                notes.append(note)
            if entry["on_base_tip"]:
                notes.append("already at the tip of `{}`".format(base_ref))
            if entry["similar_existing"]:
                similar_units = entry.get("similar_existing_units") or {}
                notes.append("similar existing: " + ", ".join(
                    "`{}`".format(s) + (" ({})".format(code_list(similar_units[s]))
                                        if similar_units.get(s) else "")
                    for s in entry["similar_existing"]))
            anchor_text = link(report, anchor)
            if differ:
                anchor_text += " (`{}`)".format(anchor.get("units", "?"))
            lines.append("| `{}` | {} | {} | {} | {} |".format(
                entry["name"], units_cell(entry),
                anchor_text,
                location_list(report, others, with_units=differ) or "",
                "; ".join(notes),
            ))
        # Drop rows from the end until the comment fits; the job summary
        # and the artifact still carry every entry.
        while len("\n".join(lines)) > MAX_MARKDOWN_CHARS and len(lines) > table_start + 1:
            lines.pop()
        hidden = len(new) - (len(lines) - table_start)
        if hidden > 0:
            lines.append("| and {} more | | | | see the workflow job summary |".format(hidden))
        lines.append("")

    if removed:
        lines.append("<details><summary>Removed standard names ({})</summary>".format(len(removed)))
        lines.append("")
        lines.append("| Removed standard name | Units | Was in | Notes |")
        lines.append("|---|---|---|---|")
        new_set = {entry["name"] for entry in new}
        for entry in removed:
            where = "<br>".join(
                "`{}:{}`".format(o["path"], o["line"]) for o in entry["occurrences"][:6])
            if len(entry["occurrences"]) > 6:
                where += "<br>and {} more".format(len(entry["occurrences"]) - 6)
            note = ""
            replacement = entry["likely_renamed_to"]
            if replacement in new_set:
                note = "likely renamed to `{}`".format(replacement)
            elif replacement:
                note = "likely replaced by the existing name `{}`".format(replacement)
            lines.append("| `{}` | {} | {} | {} |".format(
                entry["name"], units_cell(entry),
                where, note))
        lines.append("")
        lines.append("</details>")
        lines.append("")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------
def cmd_report(args: argparse.Namespace) -> int:
    report = build_report(args)
    validate_report(report)
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "report.json"), "w") as handle:
        json.dump(report, handle, indent=1)
        handle.write("\n")
    with open(os.path.join(args.out_dir, "report.md"), "w") as handle:
        handle.write(render_markdown(report))
    print("{} new and {} removed standard names (head {} vs merge base {}); "
          "report written to {}".format(
              len(report["new"]), len(report["removed"]),
              short(report["head_sha"]), short(report["merge_base_sha"]), args.out_dir))
    return 0


def cmd_render(args: argparse.Namespace) -> int:
    with open(args.report) as handle:
        report = json.load(handle)
    validate_report(report)
    text = render_markdown(report, with_marker=args.comment)
    if args.out:
        with open(args.out, "w") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)
    return 0


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    rep = sub.add_parser("report", help="compare two revisions and write the report files")
    rep.add_argument("--base", required=True, help="base branch revision, e.g. origin/main")
    rep.add_argument("--head", required=True, help="pull request head revision")
    rep.add_argument("--merge-base", help="override the computed merge base")
    rep.add_argument("--base-ref", help="base branch name for messages (default: from --base)")
    rep.add_argument("--repo", help="owner/name, enables links to the head commit")
    rep.add_argument("--pr", type=int, help="pull request number, recorded in report.json")
    rep.add_argument("--out-dir", default="stdnames_report")
    rep.set_defaults(func=cmd_report)

    ren = sub.add_parser("render", help="validate report.json and print its markdown")
    ren.add_argument("report")
    ren.add_argument("--comment", action="store_true",
                     help="prefix the hidden marker used to find the sticky comment")
    ren.add_argument("--out", help="write to this file instead of stdout")
    ren.set_defaults(func=cmd_render)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
