#!/usr/bin/env python3
"""Extract recorded seeds from SZ simulation folders into one readable text file.

Source defaults to the paper's original simulation archive; output defaults
 to the SZcerebellumModel folder on MyPassport. Reads files only; never modifies
simulation inputs or results. Requires PyYAML (included in the BSB environment).
"""
from __future__ import annotations

import argparse
import json
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

DEFAULT_SOURCE = Path(
    "/home/alberto/cerebellum/cerebellum/ref_sims/paperModelingSZ1_reproducibility"
)
DEFAULT_OUTPUT = Path(
    "/media/alberto/MyPassportAlbertoColoreROSSO/PCDESKTOPUNIPV/"
    "SZcerebellumModel/seeds_SZ.txt"
)
IGNORE_DIRS = {".git", ".ipynb_checkpoints", "__pycache__", ".venv", "BSB_venv"}
SEED_KEY = re.compile(r"(^|[_\-])(?:seed|seeds|rngseed|randomseed|grngseed)(?:$|[_\-])", re.I)
FALLBACK_LINE = re.compile(r"^\s*([\w-]*seed[\w-]*)\s*:\s*(.*?)\s*(?:#.*)?$", re.I)
RUN_NAME = re.compile(r"^20\d{12,}_")


def show(value: Any) -> str:
    if isinstance(value, (dict, list)):
        text = json.dumps(value, ensure_ascii=False, default=str)
    else:
        text = str(value)
    return text[:1500] + (" ...[truncated]" if len(text) > 1500 else "")


def find_seed_fields(obj: Any, key_path: str = ""):
    if isinstance(obj, dict):
        for key, value in obj.items():
            key_s = str(key)
            path = f"{key_path}.{key_s}" if key_path else key_s
            if SEED_KEY.search(key_s):
                yield path, value
            else:
                yield from find_seed_fields(value, path)
    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            yield from find_seed_fields(value, f"{key_path}[{i}]")


def read_yaml(path: Path):
    try:
        with path.open("r", encoding="utf-8") as f:
            return list(yaml.safe_load_all(f)), None
    except (OSError, UnicodeError, yaml.YAMLError) as error:
        return [], f"{type(error).__name__}: {error}"


def read_json(path: Path):
    try:
        with path.open("r", encoding="utf-8") as f:
            return [json.load(f)], None
    except (OSError, UnicodeError, ValueError) as error:
        return [], f"{type(error).__name__}: {error}"


def yaml_text_fallback(path: Path):
    """Best-effort lines when a YAML file cannot be parsed (not authoritative)."""
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line_no, line in enumerate(f, 1):
                match = FALLBACK_LINE.match(line)
                if match:
                    yield f"line {line_no}, {match.group(1)}", match.group(2)
    except OSError:
        return


def relevant_files(root: Path):
    for dirname, dirs, files in os.walk(root):
        dirs[:] = sorted(d for d in dirs if d not in IGNORE_DIRS)
        for filename in sorted(files):
            path = Path(dirname) / filename
            if path.suffix.lower() in {".yaml", ".yml", ".json"}:
                yield path


def deg_dirs(run: Path):
    paths = [p for p in run.iterdir() if p.is_dir() and "_deg_" in p.name]

    def order(p: Path):
        prefix = p.name.split("_deg_", 1)[0]
        return (int(prefix) if prefix.isdigit() else 10**9, p.name)

    return sorted(paths, key=order)


def scan_run(run: Path, source: Path):
    lines = ["=" * 90, f"RUN: {run.relative_to(source)}", "=" * 90]
    info = run / "info.yaml"
    recorded = []
    if info.is_file():
        docs, error = read_yaml(info)
        if error:
            lines.append(f"WARNING: cannot parse info.yaml: {error}")
        else:
            for doc in docs:
                if isinstance(doc, dict) and "seeds" in doc:
                    seeds = doc["seeds"]
                    if isinstance(seeds, list):
                        recorded = seeds
                    else:
                        lines.append(f"info.yaml :: seeds = {show(seeds)}")
    else:
        lines.append("info.yaml: not found")

    levels = deg_dirs(run)
    if recorded:
        lines.append(f"\nSEEDS IN info.yaml: {len(recorded)} (degeneration folders: {len(levels)})")
        if len(recorded) != len(levels):
            lines.append("WARNING: number of seeds != number of folders; mapping is positional only.")
        for i, seed in enumerate(recorded):
            step = levels[i].name if i < len(levels) else f"index_{i:03d} (no folder)"
            lines.append(f"  {step} -> {show(seed)}")
    else:
        lines.append("\nNo positional seed list found in info.yaml.")

    lines.append("\nSEEDS DECLARED IN YAML/JSON CONFIGURATIONS:")
    total = 0
    files_checked = 0
    parse_errors = 0
    for path in relevant_files(run):
        if path == info:
            continue  # already reported above
        files_checked += 1
        if path.stat().st_size > 20 * 1024 * 1024:
            lines.append(f"WARNING: skipped large file (>20MB): {path.relative_to(run)}")
            continue
        docs, error = read_json(path) if path.suffix.lower() == ".json" else read_yaml(path)
        if error:
            parse_errors += 1
            for field, value in yaml_text_fallback(path) if path.suffix.lower() != ".json" else []:
                lines.append(f"  {path.relative_to(run)} :: {field} = {value} [text fallback]")
                total += 1
            continue
        for i, doc in enumerate(docs):
            for field, value in find_seed_fields(doc):
                doc_label = f"document[{i}]." if len(docs) > 1 else ""
                lines.append(f"  {path.relative_to(run)} :: {doc_label}{field} = {show(value)}")
                total += 1
    if not total:
        lines.append("  None found.")
    lines.extend(["", f"Files checked: {files_checked} | seed fields found: {total} | parse errors: {parse_errors}", ""])
    return lines, len(recorded), total


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Collect seeds from SZ simulation folders without altering any data."
    )
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE,
                        help="Root containing the historical SZ run directories")
    parser.add_argument("--output", "-o", type=Path, default=DEFAULT_OUTPUT,
                        help="Destination for the text report on MyPassport")
    args = parser.parse_args()
    source = args.source.expanduser().resolve()
    output = args.output.expanduser()
    if not source.is_dir():
        parser.error(f"Source folder does not exist: {source}")

    runs = sorted(
        p for p in source.iterdir()
        if p.is_dir() and ((p / "info.yaml").is_file() or RUN_NAME.match(p.name))
    )
    if not runs:
        parser.error(f"No simulation run folders were found in: {source}")

    lines = [
        "SEED REPORT — SCHIZOPHRENIA CEREBELLUM MODEL",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Source: {source}",
        "NOTE: info.yaml seeds are recorded by the launcher, not necessarily used by NEST.",
        "NOTE: configuration fields provide independent evidence of the configured seed.",
        "",
    ]
    n_info = n_config = 0
    for run in runs:
        section, a, b = scan_run(run, source)
        lines.extend(section)
        n_info += a
        n_config += b
    lines.extend(["=" * 90,
                  f"SUMMARY: {len(runs)} runs | {n_info} launcher seeds | {n_config} configuration seed fields",
                  ""])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")
    print(f"Report saved: {output}")
    print(f"Runs: {len(runs)}, launcher seeds: {n_info}, configuration fields: {n_config}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
