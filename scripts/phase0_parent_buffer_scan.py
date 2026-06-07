"""
Phase 0 parent-buffer scan for the Gaia DR3 slow-Vgrf catalogue.

This script scans the local Gaia DR3 CSV mirror and writes a broad legacy
low-Vgrf buffer from the full Gaia 6D parent sample. The buffer is intentionally
wider than the final 25 km/s catalogue threshold so that later phases can prove
that the final probabilistic catalogue was not limited by the old <25 km/s
preselection.

The scan is designed for an external mechanical drive:
  - checkpoint/resume
  - short sleeps between files
  - longer cooling pauses every batch
  - compact per-file progress and a JSON manifest

Default output:
  private_outputs/parent_scan/gaia_parent_buffer_vgrf200_full.csv
  private_outputs/parent_scan/gaia_parent_buffer_scan_state.json
  private_outputs/parent_scan/gaia_parent_buffer_scan_manifest.json
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


BUNDLE = Path(__file__).resolve().parents[1]
REPO = BUNDLE.parents[1]
DEFAULT_DATA_DIR = BUNDLE / "private_inputs" / "gaia_dr3_source_mirror"
DEFAULT_OUT_DIR = BUNDLE / "private_outputs" / "parent_scan"

K = 4.74047
U_SUN = 11.1
V_SUN = 12.24
W_SUN = 7.25
V_CIRC = 229.0  # WP1/O1 (Phase-2) adopt 229; early screen not re-run

AG = np.array(
    [
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [+0.4941094279, -0.4448296300, +0.7469822445],
        [-0.8676661490, -0.1980763734, +0.4559837762],
    ]
)

COUNT_BINS = [25, 50, 75, 100, 125, 150, 175, 200]


def utc_now() -> str:
    return datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def count_header_lines(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                count += 1
            else:
                break
    return count


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    if len(df.columns) and df.columns[0].startswith("\ufeff"):
        df.columns = [df.columns[0].replace("\ufeff", "")] + list(df.columns[1:])
    return df


def compute_legacy_vgrf(df: pd.DataFrame, min_parallax_over_error: float) -> pd.DataFrame:
    required = ["source_id", "ra", "dec", "parallax", "pmra", "pmdec", "radial_velocity"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"missing required columns: {missing}")

    mask = (
        df["radial_velocity"].notna()
        & df["parallax"].notna()
        & df["pmra"].notna()
        & df["pmdec"].notna()
        & (df["parallax"] > 0)
    )
    if "parallax_over_error" in df.columns:
        mask &= df["parallax_over_error"] >= min_parallax_over_error
    elif "parallax_error" in df.columns:
        mask &= (df["parallax"] / df["parallax_error"]) >= min_parallax_over_error

    usable = df.loc[mask].copy()
    if usable.empty:
        return usable

    d_kpc = 1.0 / usable["parallax"].to_numpy(dtype=float)
    ra_rad = np.radians(usable["ra"].to_numpy(dtype=float))
    dec_rad = np.radians(usable["dec"].to_numpy(dtype=float))
    pmra = usable["pmra"].to_numpy(dtype=float)
    pmdec = usable["pmdec"].to_numpy(dtype=float)
    rv = usable["radial_velocity"].to_numpy(dtype=float)

    v_ra = K * pmra * d_kpc
    v_dec = K * pmdec * d_kpc

    cos_ra = np.cos(ra_rad)
    sin_ra = np.sin(ra_rad)
    cos_dec = np.cos(dec_rad)
    sin_dec = np.sin(dec_rad)

    vx_eq = rv * cos_ra * cos_dec - v_ra * sin_ra - v_dec * cos_ra * sin_dec
    vy_eq = rv * sin_ra * cos_dec + v_ra * cos_ra - v_dec * sin_ra * sin_dec
    vz_eq = rv * sin_dec + v_dec * cos_dec

    v_gal = AG @ np.vstack([vx_eq, vy_eq, vz_eq])
    u_lsr = v_gal[0] + U_SUN
    v_lsr = v_gal[1] + V_SUN
    w_lsr = v_gal[2] + W_SUN
    v_grf = v_lsr + V_CIRC
    v_total_grf = np.sqrt(u_lsr**2 + v_grf**2 + w_lsr**2)

    usable["legacy_distance_pc"] = d_kpc * 1000.0
    usable["legacy_U_lsr"] = u_lsr
    usable["legacy_V_lsr"] = v_lsr
    usable["legacy_W_lsr"] = w_lsr
    usable["legacy_V_grf"] = v_grf
    usable["legacy_v_total_grf"] = v_total_grf
    return usable


def load_state(path: Path) -> dict[str, Any]:
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    return {
        "created_utc": utc_now(),
        "completed_files": [],
        "total_rows": 0,
        "total_parent_6d": 0,
        "total_buffer_rows": 0,
        "bin_counts": {f"lt_{bin_value}": 0 for bin_value in COUNT_BINS},
        "per_file": [],
        "done": False,
    }


def save_json(path: Path, data: dict[str, Any]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
    tmp.replace(path)


def append_csv(path: Path, rows: pd.DataFrame) -> None:
    write_header = not path.exists()
    rows.to_csv(path, mode="a", index=False, header=write_header)


def run_scan(args: argparse.Namespace) -> dict[str, Any]:
    data_dir = Path(args.data_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    output_csv = out_dir / f"gaia_parent_buffer_vgrf{int(args.buffer_kms)}_full.csv"
    state_path = out_dir / "gaia_parent_buffer_scan_state.json"
    manifest_path = out_dir / "gaia_parent_buffer_scan_manifest.json"
    log_path = out_dir / "gaia_parent_buffer_scan.log"

    files = sorted(Path(p) for p in glob.glob(str(data_dir / "GaiaSource_*.csv")))
    if args.max_files is not None:
        files = files[: args.max_files]
    if not files:
        raise FileNotFoundError(f"No GaiaSource_*.csv files found under {data_dir}")

    state = load_state(state_path)
    if state.get("done") and not args.force:
        print(f"Scan is already marked done in {state_path}")
        return state

    completed = set(state.get("completed_files", []))
    started = time.time()

    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"[{utc_now()}] starting/resuming scan over {len(files)} files\n")
        log.write(f"[{utc_now()}] buffer_kms={args.buffer_kms}, data_dir={data_dir}\n")

    for file_index, path in enumerate(files, start=1):
        name = path.name
        if name in completed:
            continue

        file_started = time.time()
        try:
            skiprows = count_header_lines(path)
            n_rows = 0
            n_parent = 0
            n_buffer = 0
            counts = {f"lt_{bin_value}": 0 for bin_value in COUNT_BINS}

            reader = pd.read_csv(
                path,
                skiprows=skiprows,
                low_memory=False,
                chunksize=args.read_chunksize,
            )
            for df in reader:
                df = clean_columns(df)
                velocity_df = compute_legacy_vgrf(df, args.min_parallax_over_error)

                if velocity_df.empty:
                    buffer_rows = velocity_df
                    v = np.array([], dtype=float)
                else:
                    buffer_rows = velocity_df.loc[
                        velocity_df["legacy_v_total_grf"] < args.buffer_kms
                    ].copy()
                    if not buffer_rows.empty:
                        buffer_rows.insert(0, "parent_scan_file", name)
                        append_csv(output_csv, buffer_rows)
                    v = velocity_df["legacy_v_total_grf"].to_numpy(dtype=float)

                n_rows += int(len(df))
                n_parent += int(len(velocity_df))
                n_buffer += int(len(buffer_rows))
                for bin_value in COUNT_BINS:
                    counts[f"lt_{bin_value}"] += int(np.sum(v < bin_value))

            state["total_rows"] += n_rows
            state["total_parent_6d"] += n_parent
            state["total_buffer_rows"] += n_buffer
            for key, value in counts.items():
                state["bin_counts"][key] = int(state["bin_counts"].get(key, 0) + value)

            elapsed = time.time() - file_started
            state["per_file"].append(
                {
                    "file": name,
                    "rows": int(n_rows),
                    "parent_6d": int(n_parent),
                    "buffer_rows": int(n_buffer),
                    "bin_counts": counts,
                    "seconds": round(elapsed, 3),
                }
            )
            completed.add(name)
            state["completed_files"] = sorted(completed)
            state["updated_utc"] = utc_now()
            state["done"] = False

            if len(completed) % args.checkpoint_every == 0:
                save_json(state_path, state)

            done = len(completed)
            wall = time.time() - started
            rate = done / max(wall, 1.0)
            remaining = (len(files) - done) / max(rate, 0.001)
            msg = (
                f"[{done}/{len(files)}] {done * 100 / len(files):5.1f}% "
                f"{name} rows={n_rows:,} parent6d={n_parent:,} "
                f"buffer<{args.buffer_kms:g}={n_buffer:,} "
                f"total_buffer={state['total_buffer_rows']:,} "
                f"{elapsed:.1f}s ETA={timedelta(seconds=int(remaining))}"
            )
            print(msg, flush=True)
            with log_path.open("a", encoding="utf-8") as log:
                log.write(f"[{utc_now()}] {msg}\n")

        except Exception as exc:
            err = f"ERROR {name}: {type(exc).__name__}: {exc}"
            print(err, flush=True)
            with log_path.open("a", encoding="utf-8") as log:
                log.write(f"[{utc_now()}] {err}\n")

        save_json(state_path, state)

        if args.throttle_seconds > 0:
            time.sleep(args.throttle_seconds)
        if (
            args.batch_pause_every > 0
            and len(completed) > 0
            and len(completed) % args.batch_pause_every == 0
            and args.batch_pause_seconds > 0
        ):
            print(
                f"Cooling pause: {args.batch_pause_seconds:g}s "
                f"after {len(completed)} completed files",
                flush=True,
            )
            time.sleep(args.batch_pause_seconds)

    state["done"] = True
    state["completed_utc"] = utc_now()
    state["output_csv"] = str(output_csv)
    state["manifest_path"] = str(manifest_path)
    state["parameters"] = vars(args)
    state["files_seen"] = len(files)
    save_json(state_path, state)
    save_json(manifest_path, state)
    return state


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default=str(DEFAULT_DATA_DIR))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--buffer-kms", type=float, default=200.0)
    parser.add_argument("--min-parallax-over-error", type=float, default=5.0)
    parser.add_argument("--throttle-seconds", type=float, default=0.5)
    parser.add_argument("--batch-pause-every", type=int, default=50)
    parser.add_argument("--batch-pause-seconds", type=float, default=5.0)
    parser.add_argument("--checkpoint-every", type=int, default=5)
    parser.add_argument("--read-chunksize", type=int, default=100_000)
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    state = run_scan(args)
    print(
        "Scan summary: "
        f"files={len(state.get('completed_files', []))}, "
        f"rows={state.get('total_rows', 0):,}, "
        f"parent6d={state.get('total_parent_6d', 0):,}, "
        f"buffer={state.get('total_buffer_rows', 0):,}, "
        f"bins={state.get('bin_counts', {})}"
    )


if __name__ == "__main__":
    main()
