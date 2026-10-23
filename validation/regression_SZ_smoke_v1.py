#!/usr/bin/env python3
"""Numerical regression checker for SZ BASAL / EVOKED / AUTOR smoke CSVs.

Uses only the Python standard library. Never executes notebooks, simulations,
validator scripts, or pickle files. Baselines are provisional smoke-test outputs,
not independently established scientific ground truth.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DATASETS = {
    "BASAL": ("basal_delta_firing_smoke.csv", "basal_bridge_check.json"),
    "EVOKED": ("evoked_smoke_metrics.csv", "evoked_bridge_check.json"),
    "AUTOR": ("autor_smoke_delta_firing.csv", "autor_smoke_check.json"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.reader(file)
        try:
            fields = next(reader)
        except StopIteration as error:
            raise ValueError(f"Empty CSV: {path}") from error
        if not fields or any(not field for field in fields):
            raise ValueError(f"Empty column name: {path}")
        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate column names: {path}")
        records = []
        for lineno, record in enumerate(reader, start=2):
            if len(record) != len(fields):
                raise ValueError(
                    f"CSV {path}, line {lineno}: expected {len(fields)} cells, "
                    f"found {len(record)}"
                )
            records.append(record)
    if not records:
        raise ValueError(f"CSV has a header but no data: {path}")
    return fields, records


def check_json(path: Path) -> None:
    with path.open("r", encoding="utf-8") as file:
        json.load(file)


def git_revision() -> str | None:
    try:
        r = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True, timeout=3
        )
        return r.stdout.strip()
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None


def input_files(args: argparse.Namespace) -> dict[str, Path]:
    return {
        "BASAL": args.basal.expanduser().resolve(),
        "EVOKED": args.evoked.expanduser().resolve(),
        "AUTOR": args.autor.expanduser().resolve(),
    }


def capture(args: argparse.Namespace) -> int:
    sources = input_files(args)
    baseline = args.baseline_dir.expanduser().resolve()
    for name, path in sources.items():
        if not path.is_file():
            raise FileNotFoundError(f"{name}: input CSV missing: {path}")
        load_csv(path)

    # Refuse overwrite, even if a prior baseline directory is empty.
    if baseline.exists():
        raise FileExistsError(
            f"Baseline already exists: {baseline}. "
            "Use another versioned directory; don't silently overwrite it."
        )

    # Avoid creating a partial baseline when checks cannot be parsed.
    for name, path in sources.items():
        _, check_name = DATASETS[name]
        check_path = path.parent / check_name
        if check_path.is_file():
            check_json(check_path)

    baseline.mkdir(parents=True, exist_ok=False)
    info: dict[str, Any] = {
        "format_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "git_revision": git_revision(),
        "description": (
            "Provisional numerical baselines captured from already PASS smoke outputs; "
            "not independent evidence of agreement with original historical publications."
        ),
        "datasets": {},
    }

    for name, original in sources.items():
        filename, check_name = DATASETS[name]
        fields, rows = load_csv(original)
        snap = baseline / filename
        shutil.copy2(original, snap)
        check_path = original.parent / check_name
        check_sha = None
        if check_path.is_file():
            shutil.copy2(check_path, baseline / check_name)
            check_sha = sha256(baseline / check_name)
        info["datasets"][name] = {
            "csv_filename": filename,
            "csv_sha256": sha256(snap),
            "check_json_filename": check_name if check_sha else None,
            "check_json_sha256": check_sha,
            "source_csv": str(original),
            "source_check_json": str(check_path) if check_sha else None,
            "row_count": len(rows),
            "columns": fields,
        }
        print(f"CAPTURE {name}: {len(rows)} rows × {len(fields)} columns -> {snap}")

    manifest = baseline / "manifest.json"
    manifest.write_text(json.dumps(info, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"BASELINE CAPTURED: {manifest}")
    print("CAUTION: Smoke-test baseline, not externally validated scientific ground truth.")
    return 0


def num_value(s: str) -> float | None:
    if not s or not s.strip():
        return None
    try:
        return float(s)
    except ValueError:
        return None


def compare_cell(expected: str, observed: str, rtol: float, atol: float) -> tuple[bool, float | None]:
    x, y = num_value(expected), num_value(observed)
    if x is not None and y is not None:
        if math.isnan(x) or math.isnan(y):
            return math.isnan(x) and math.isnan(y), None
        if math.isinf(x) or math.isinf(y):
            return x == y, None
        delta = abs(y - x)
        return math.isclose(x, y, rel_tol=rtol, abs_tol=atol), delta
    return expected == observed, None


def compare_data(
    expected_file: Path, observed_file: Path, rtol: float, atol: float
) -> dict[str, Any]:
    h0, d0 = load_csv(expected_file)
    h1, d1 = load_csv(observed_file)
    result: dict[str, Any] = {
        "passed": True,
        "baseline_rows": len(d0),
        "observed_rows": len(d1),
        "baseline_columns": h0,
        "observed_columns": h1,
        "cells_compared": 0,
        "numeric_cells_compared": 0,
        "max_abs_difference": 0.0,
        "mismatch_count": 0,
        "mismatch_examples": [],
    }
    if h0 != h1:
        result["passed"] = False
        result["mismatch_count"] += 1
        result["mismatch_examples"].append({"issue": "Column names/order differ"})
        return result
    if len(d0) != len(d1):
        result["passed"] = False
        result["mismatch_count"] += 1
        result["mismatch_examples"].append({"issue": "Row count differs"})

    for idx, (row0, row1) in enumerate(zip(d0, d1), start=2):
        for col, v0, v1 in zip(h0, row0, row1):
            result["cells_compared"] += 1
            x, y = num_value(v0), num_value(v1)
            if x is not None and y is not None:
                result["numeric_cells_compared"] += 1
            same, difference = compare_cell(v0, v1, rtol, atol)
            if difference is not None:
                result["max_abs_difference"] = max(result["max_abs_difference"], difference)
            if not same:
                result["passed"] = False
                result["mismatch_count"] += 1
                if len(result["mismatch_examples"]) < 12:
                    result["mismatch_examples"].append({
                        "csv_line": idx,
                        "column": col,
                        "baseline": v0,
                        "observed": v1,
                        "abs_difference": difference,
                    })
    return result


def compare(args: argparse.Namespace) -> int:
    baseline = args.baseline_dir.expanduser().resolve()
    manifest_path = baseline / "manifest.json"
    with manifest_path.open("r", encoding="utf-8") as file:
        info = json.load(file)
    if info.get("format_version") != 1:
        raise ValueError("Unsupported baseline format_version")

    report: dict[str, Any] = {
        "format_version": 1,
        "checked_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "baseline_manifest": str(manifest_path),
        "rtol": args.rtol,
        "atol": args.atol,
        "row_matching": "strict file order (header and row positions)",
        "datasets": {},
    }
    inputs = input_files(args)
    passed = True
    for name, observed in inputs.items():
        record = info["datasets"].get(name)
        if not isinstance(record, dict):
            raise ValueError(f"Manifest lacks dataset: {name}")
        filename, check_name = DATASETS[name]
        saved = baseline / filename
        if sha256(saved) != record.get("csv_sha256"):
            raise ValueError(f"{name}: baseline CSV checksum changed; do not compare")
        if record.get("check_json_sha256"):
            saved_check = baseline / check_name
            if sha256(saved_check) != record["check_json_sha256"]:
                raise ValueError(f"{name}: baseline JSON checksum changed; do not compare")
        if not observed.is_file():
            raise FileNotFoundError(f"{name}: fresh CSV missing: {observed}")
        if str(observed) == record.get("source_csv") or observed == saved:
            raise ValueError(
                f"{name}: comparison uses the captured source or baseline CSV itself. "
                "Run the validator again in a distinct output directory first."
            )
        check = observed.parent / check_name
        if not check.is_file():
            raise FileNotFoundError(f"{name}: new diagnostic JSON missing: {check}")
        check_json(check)
        stats = compare_data(saved, observed, args.rtol, args.atol)
        stats["observed_csv"] = str(observed)
        stats["observed_check_json"] = str(check)
        report["datasets"][name] = stats
        if stats["passed"]:
            print(
                f"PASS {name}: {stats['observed_rows']} rows, "
                f"{len(stats['observed_columns'])} columns, "
                f"{stats['numeric_cells_compared']} numeric cells; "
                f"max |Δ| = {stats['max_abs_difference']:.8g}"
            )
        else:
            passed = False
            print(f"FAIL {name}: {stats['mismatch_count']} differences")
            for example in stats["mismatch_examples"]:
                print(f"  {example}")

    report["passed"] = passed
    if args.report:
        dest = args.report.expanduser().resolve()
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"REPORT: {dest}")
    print(f"SZ NUMERICAL REGRESSION: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("capture", "compare"):
        cmd = commands.add_parser(name)
        cmd.add_argument("--baseline-dir", type=Path, required=True)
        cmd.add_argument("--basal", type=Path, required=True, help="Path to BASAL CSV")
        cmd.add_argument("--evoked", type=Path, required=True, help="Path to EVOKED CSV")
        cmd.add_argument("--autor", type=Path, required=True, help="Path to AUTOR CSV")
        if name == "compare":
            cmd.add_argument("--rtol", type=float, default=1e-6)
            cmd.add_argument("--atol", type=float, default=1e-8)
            cmd.add_argument("--report", type=Path)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "compare":
        if not (math.isfinite(args.rtol) and math.isfinite(args.atol)):
            parser.error("rtol and atol must be finite")
        if args.rtol < 0 or args.atol < 0:
            parser.error("rtol and atol must be nonnegative")
    try:
        return capture(args) if args.command == "capture" else compare(args)
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
