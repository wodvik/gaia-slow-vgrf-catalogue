"""Phase 14X -- full expanded-catalogue Monte Carlo orbit propagation.

This is the expanded replacement for the older Phase 6A orbit-MC pass.
It runs a high-sample static-potential orbit Monte Carlo for every
Tier A+B+C star in the expanded catalogue.

Default plan:
  - 5,000 realisations for each expanded Tier A+B+C star.
  - 10,000-realisation convergence rerun for a fixed random subset.
  - Chunked output written to disk after every chunk.

Run under WSL, where AGAMA is installed:

    python3 scripts/compute_mc_orbits.py
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import time
from pathlib import Path

import agama
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table


if platform.system().lower() == "linux":
    REPO = Path(__file__).resolve().parents[1]
    DEFAULT_INPUT = Path("/mnt/d/GAIA/parent_scan/expanded_candidates_mc_tiered.csv")
else:
    REPO = Path(__file__).resolve().parents[1]
    DEFAULT_INPUT = Path("D:/GAIA/parent_scan/expanded_candidates_mc_tiered.csv")

CONFIG = yaml.safe_load((REPO / "config.yml").read_text())
OUT = REPO / "analysis_products/expanded_orbit_mc"
WORK = REPO / "phase3_agama/_hunter24_workdir"
SEED = int(CONFIG["mc"]["random_seed"])
GYR = 1.0 / 0.9778

agama.setUnits(length=1, mass=1, velocity=1)


def log(message: str) -> None:
    print(f"[14X-expanded-orbit-MC t={time.time() - T0:8.1f}s] {message}", flush=True)


def galcen_frame() -> coord.Galactocentric:
    sv = CONFIG["solar_variants"]["default"]
    return coord.Galactocentric(
        galcen_distance=sv["R0_kpc"] * u.kpc,
        z_sun=sv["z_sun_pc"] * u.pc,
        galcen_v_sun=coord.CartesianDifferential(
            sv["U_kms"] * u.km / u.s,
            (sv["Vc_kms"] + sv["V_kms"]) * u.km / u.s,
            sv["W_kms"] * u.km / u.s,
        ),
    )


def cov3(df: pd.DataFrame) -> np.ndarray:
    sig = np.stack(
        [
            df["parallax_error"].to_numpy(dtype=float),
            df["pmra_error"].to_numpy(dtype=float),
            df["pmdec_error"].to_numpy(dtype=float),
        ],
        axis=-1,
    )
    rho = np.zeros((len(df), 3, 3), dtype=float)
    for i in (0, 1, 2):
        rho[:, i, i] = 1.0
    rho[:, 0, 1] = rho[:, 1, 0] = df["parallax_pmra_corr"].to_numpy(dtype=float)
    rho[:, 0, 2] = rho[:, 2, 0] = df["parallax_pmdec_corr"].to_numpy(dtype=float)
    rho[:, 1, 2] = rho[:, 2, 1] = df["pmra_pmdec_corr"].to_numpy(dtype=float)
    return rho * (sig[:, :, None] * sig[:, None, :])


def chol3(cov: np.ndarray) -> np.ndarray:
    out = np.zeros_like(cov)
    eye = np.eye(3)
    for i in range(cov.shape[0]):
        try:
            out[i] = np.linalg.cholesky(cov[i])
        except np.linalg.LinAlgError:
            out[i] = np.linalg.cholesky(cov[i] + 1e-9 * eye)
    return out


def split_normal(med: np.ndarray, lo: np.ndarray, hi: np.ndarray, n_samp: int, rng: np.random.Generator) -> np.ndarray:
    sig_lo = np.maximum(med - lo, 1.0)
    sig_hi = np.maximum(hi - med, 1.0)
    z = rng.standard_normal((len(med), n_samp))
    sigma = np.where(z < 0.0, sig_lo[:, None], sig_hi[:, None])
    return np.maximum(med[:, None] + z * sigma, 1.0)


def sample_ics(df: pd.DataFrame, n_samp: int, rng: np.random.Generator) -> np.ndarray:
    mu = np.stack(
        [
            df["parallax_zpcorr"].to_numpy(dtype=float),
            df["pmra"].to_numpy(dtype=float),
            df["pmdec"].to_numpy(dtype=float),
        ],
        axis=-1,
    )
    z = rng.standard_normal((len(df), n_samp, 3))
    ast = mu[:, None, :] + np.einsum("nij,nsj->nsi", chol3(cov3(df)), z)
    rv = (
        df["radial_velocity"].to_numpy(dtype=float)[:, None]
        + rng.standard_normal((len(df), n_samp))
        * df["radial_velocity_error"].to_numpy(dtype=float)[:, None]
    )
    dist = split_normal(
        df["dist_pc_final_screen"].to_numpy(dtype=float),
        df["dist_lo_pc_final_screen"].to_numpy(dtype=float),
        df["dist_hi_pc_final_screen"].to_numpy(dtype=float),
        n_samp,
        rng,
    )
    ra = np.broadcast_to(df["ra"].to_numpy(dtype=float)[:, None], (len(df), n_samp))
    dec = np.broadcast_to(df["dec"].to_numpy(dtype=float)[:, None], (len(df), n_samp))

    icrs = coord.SkyCoord(
        ra=ra.ravel() * u.deg,
        dec=dec.ravel() * u.deg,
        distance=dist.ravel() * u.pc,
        pm_ra_cosdec=ast[:, :, 1].ravel() * u.mas / u.yr,
        pm_dec=ast[:, :, 2].ravel() * u.mas / u.yr,
        radial_velocity=rv.ravel() * u.km / u.s,
        frame="icrs",
    )
    g = icrs.transform_to(galcen_frame())
    ic = np.column_stack(
        [
            g.x.to_value(u.kpc),
            g.y.to_value(u.kpc),
            g.z.to_value(u.kpc),
            g.v_x.to_value(u.km / u.s),
            g.v_y.to_value(u.km / u.s),
            g.v_z.to_value(u.km / u.s),
        ]
    )
    return ic.reshape(len(df), n_samp, 6)


def integrate_orbits(ic_flat: np.ndarray, pot, n_steps: int) -> tuple[np.ndarray, ...]:
    res = agama.orbit(potential=pot, ic=ic_flat, time=4.0 * GYR, trajsize=n_steps)
    n = len(ic_flat)
    rperi = np.empty(n)
    rapo = np.empty(n)
    zmax = np.empty(n)
    ecc = np.empty(n)
    rsph_min = np.empty(n)
    for i in range(n):
        traj = np.asarray(res[i, 1])
        rcyl = np.hypot(traj[:, 0], traj[:, 1])
        rsph = np.sqrt(rcyl * rcyl + traj[:, 2] * traj[:, 2])
        rperi[i] = rcyl.min()
        rapo[i] = rcyl.max()
        zmax[i] = np.abs(traj[:, 2]).max()
        ecc[i] = (rapo[i] - rperi[i]) / (rapo[i] + rperi[i])
        rsph_min[i] = rsph.min()
    return rperi, rapo, zmax, ecc, rsph_min


def percentiles(values: np.ndarray) -> np.ndarray:
    return np.percentile(values, [16, 50, 84], axis=1).T


def summarise_chunk(df: pd.DataFrame, ic_3d: np.ndarray, pot, af, n_steps: int, n_samp: int) -> pd.DataFrame:
    flat = ic_3d.reshape(-1, 6)
    rperi, rapo, zmax, ecc, rsph = integrate_orbits(flat, pot, n_steps)
    actions = af(flat)
    n = len(df)
    arrays = {
        "R_peri_kpc": rperi.reshape(n, n_samp),
        "R_apo_kpc": rapo.reshape(n, n_samp),
        "z_max_kpc": zmax.reshape(n, n_samp),
        "ecc": ecc.reshape(n, n_samp),
        "min_r_sph_kpc": rsph.reshape(n, n_samp),
        "J_R": actions[:, 0].reshape(n, n_samp),
        "J_z": actions[:, 1].reshape(n, n_samp),
        "J_phi": actions[:, 2].reshape(n, n_samp),
    }
    out = df[["source_id", "tier", "P_vgrf_below_25", "vgrf_default_exact", "mc_realisations"]].copy()
    out["orbit_mc_realisations"] = n_samp
    for name, arr in arrays.items():
        p = percentiles(arr)
        out[f"{name}_p16"] = p[:, 0]
        out[f"{name}_p50"] = p[:, 1]
        out[f"{name}_p84"] = p[:, 2]
    return out


def load_expanded(input_csv: Path) -> pd.DataFrame:
    cols = [
        "source_id", "tier", "P_vgrf_below_25", "mc_realisations",
        "vgrf_default_exact", "ra", "dec", "parallax_zpcorr",
        "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error",
        "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr",
        "radial_velocity", "radial_velocity_error",
        "dist_pc_final_screen", "dist_lo_pc_final_screen",
        "dist_hi_pc_final_screen",
    ]
    df = pd.read_csv(input_csv, usecols=cols)
    df["tier"] = df["tier"].astype(str).str.strip()
    df["source_id"] = df["source_id"].astype("int64")
    df = df[df["tier"].isin(["A", "B", "C"])].copy()
    return df.sort_values(["tier", "source_id"]).reset_index(drop=True)


def combine_chunks(chunks_dir: Path, output_dir: Path) -> pd.DataFrame:
    parts = [pd.read_csv(path) for path in sorted(chunks_dir.glob("chunk_*.csv"))]
    if not parts:
        raise RuntimeError(f"No chunk outputs found under {chunks_dir}")
    out = pd.concat(parts, ignore_index=True).sort_values("source_id").reset_index(drop=True)
    out_csv = output_dir / "expanded_catalogue_mc_orbits.csv"
    out_fits = output_dir / "expanded_catalogue_mc_orbits.fits"
    out.to_csv(out_csv, index=False)
    Table.from_pandas(out).write(out_fits, overwrite=True)
    return out


def stat(values: pd.Series) -> dict[str, float | int]:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    return {
        "n": int(len(arr)),
        "p16": float(np.percentile(arr, 16)),
        "p50": float(np.percentile(arr, 50)),
        "p84": float(np.percentile(arr, 84)),
        "min": float(arr.min()),
        "max": float(arr.max()),
    }


def run_convergence(df: pd.DataFrame, pot, af, args: argparse.Namespace) -> dict:
    conv_dir = OUT / "convergence_chunks"
    conv_dir.mkdir(parents=True, exist_ok=True)
    rng_pick = np.random.default_rng(SEED + 14900)
    n_test = min(args.convergence_stars, len(df))
    test_idx_path = OUT / "convergence_source_ids.json"
    if test_idx_path.exists() and not args.force_convergence:
        source_ids = json.loads(test_idx_path.read_text())["source_ids"]
        test = df[df["source_id"].isin(source_ids)].copy().reset_index(drop=True)
    else:
        source_ids = df.iloc[rng_pick.choice(len(df), size=n_test, replace=False)]["source_id"].astype(int).tolist()
        test_idx_path.write_text(json.dumps({"source_ids": source_ids}, indent=2))
        test = df[df["source_id"].isin(source_ids)].copy().reset_index(drop=True)

    log(f"convergence: {len(test)} stars x {args.n_convergence:,} realisations")
    rng = np.random.default_rng(SEED + 14901)
    for start in range(0, len(test), args.convergence_chunk_size):
        stop = min(start + args.convergence_chunk_size, len(test))
        path = conv_dir / f"conv_{start:04d}_{stop:04d}.csv"
        if path.exists() and not args.force_convergence:
            log(f"convergence chunk {start}-{stop} exists; skipping")
            continue
        log(f"convergence chunk {start}-{stop}")
        chunk = test.iloc[start:stop].reset_index(drop=True)
        ic = sample_ics(chunk, args.n_convergence, rng)
        out = summarise_chunk(chunk, ic, pot, af, args.trajsize, args.n_convergence)
        out.to_csv(path, index=False)
        log(f"wrote {path}")

    conv = pd.concat([pd.read_csv(path) for path in sorted(conv_dir.glob("conv_*.csv"))], ignore_index=True)
    conv_csv = OUT / "expanded_catalogue_mc_orbits_convergence_10000.csv"
    conv.to_csv(conv_csv, index=False)
    return {
        "n_convergence_stars": int(len(conv)),
        "n_convergence_realisations": int(args.n_convergence),
        "convergence_csv": str(conv_csv),
        "R_peri_kpc_p50": stat(conv["R_peri_kpc_p50"]),
        "ecc_p50": stat(conv["ecc_p50"]),
    }


def main() -> int:
    global OUT
    global T0
    T0 = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--n-samp", type=int, default=5000)
    parser.add_argument("--chunk-size", type=int, default=25)
    parser.add_argument("--trajsize", type=int, default=1001)
    parser.add_argument("--n-convergence", type=int, default=10000)
    parser.add_argument("--convergence-stars", type=int, default=100)
    parser.add_argument("--convergence-chunk-size", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--force-convergence", action="store_true")
    args = parser.parse_args()

    OUT = args.out_dir
    OUT.mkdir(parents=True, exist_ok=True)
    chunks_dir = OUT / "chunks"
    chunks_dir.mkdir(parents=True, exist_ok=True)

    log(f"agama {agama.__version__}; os.cpu_count={os.cpu_count()}; OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}")
    log(f"input={args.input_csv}")
    df = load_expanded(args.input_csv)
    log(f"expanded Tier A+B+C rows: {len(df)}")

    pot = agama.Potential(file=str(WORK / "MWPotentialHunter24_axi.ini"))
    af = agama.ActionFinder(pot, interp=True)

    rng = np.random.default_rng(SEED + 14000)
    for start in range(0, len(df), args.chunk_size):
        stop = min(start + args.chunk_size, len(df))
        path = chunks_dir / f"chunk_{start:04d}_{stop:04d}.csv"
        if path.exists() and not args.force:
            log(f"chunk {start}-{stop} exists; skipping")
            continue
        chunk = df.iloc[start:stop].reset_index(drop=True)
        log(f"chunk {start}-{stop}: {len(chunk)} stars x {args.n_samp:,} realisations")
        ic = sample_ics(chunk, args.n_samp, rng)
        out = summarise_chunk(chunk, ic, pot, af, args.trajsize, args.n_samp)
        out.to_csv(path, index=False)
        log(f"wrote {path}")

    full = combine_chunks(chunks_dir, OUT)
    convergence = run_convergence(df, pot, af, args)
    summary = {
        "phase": "14X",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_csv": str(args.input_csv),
        "n_tierABC": int(len(df)),
        "n_realisations_per_star": int(args.n_samp),
        "trajsize": int(args.trajsize),
        "chunk_size": int(args.chunk_size),
        "threads": {
            "os_cpu_count": os.cpu_count(),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
        },
        "outputs": {
            "fits": str(OUT / "expanded_catalogue_mc_orbits.fits"),
            "csv": str(OUT / "expanded_catalogue_mc_orbits.csv"),
            "chunks_dir": str(chunks_dir),
        },
        "headline": {
            "R_peri_p50_pc": stat(full["R_peri_kpc_p50"] * 1000.0),
            "R_apo_p50_kpc": stat(full["R_apo_kpc_p50"]),
            "ecc_p50": stat(full["ecc_p50"]),
            "z_max_p50_kpc": stat(full["z_max_kpc_p50"]),
        },
        "per_tier": {},
        "convergence": convergence,
    }
    for tier in ["A", "B", "C"]:
        sub = full[full["tier"] == tier]
        summary["per_tier"][tier] = {
            "n": int(len(sub)),
            "R_peri_p50_pc": stat(sub["R_peri_kpc_p50"] * 1000.0),
            "R_apo_p50_kpc": stat(sub["R_apo_kpc_p50"]),
            "ecc_p50": stat(sub["ecc_p50"]),
        }
    summary_path = OUT / "expanded_catalogue_mc_orbits_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    log(f"DONE; wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
