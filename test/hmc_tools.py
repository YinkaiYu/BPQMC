from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, replace
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = ROOT / "src" / "BPQMC.out"
ACCEPT_RE = re.compile(r"^\s*([A-Za-z0-9_ ]+?)\s*:\s*([0-9Ee+\-.]+)\s*$")
N_AUX = 2


@dataclass(frozen=True)
class RunConfig:
    name: str
    lattice_type: str
    rt: float
    ru1: float
    ru2: float
    nbos: int
    nlx: int
    nly: int
    ltrot: int
    beta: float
    nlx_therm: int | None = None
    nly_therm: int | None = None
    ltrot_therm: int | None = None
    nwrap: int = 2
    nbin: int = 16
    nsweep: int = 4
    shift_loc: float = 1.0
    is_tau: bool = False
    nthermal: int = 0
    is_warm: bool = False
    nwarm: int = 0
    shift_warm_1: float = 1.0
    shift_warm_2: float = 1.0
    ini_type: int = 2
    ini_ampl: float = 0.1
    ini_bias_1: float = 0.0
    ini_bias_2: float = 0.0
    ini_ham: int = 0
    ini_twist: float = 1.0e-4
    imbalance: float = 0.0

    def with_sampling(self, *, nbin: int | None = None, nsweep: int | None = None, is_warm: bool | None = None, nwarm: int | None = None) -> "RunConfig":
        return replace(
            self,
            nbin=self.nbin if nbin is None else nbin,
            nsweep=self.nsweep if nsweep is None else nsweep,
            is_warm=self.is_warm if is_warm is None else is_warm,
            nwarm=self.nwarm if nwarm is None else nwarm,
        )

    def therm_dims(self) -> tuple[int, int, int]:
        return (
            self.nlx if self.nlx_therm is None else self.nlx_therm,
            self.nly if self.nly_therm is None else self.nly_therm,
            self.ltrot if self.ltrot_therm is None else self.ltrot_therm,
        )


def default_parameter_sets() -> dict[str, RunConfig]:
    return {
        "weak_u2": RunConfig(
            name="weak_u2",
            lattice_type="kagome",
            rt=1.0,
            ru1=0.0,
            ru2=1.0,
            nbos=6,
            nlx=2,
            nly=2,
            ltrot=8,
            beta=1.6,
            nwrap=2,
            nbin=400,
            nsweep=8,
            nthermal=200,
        ),
        "mixed_u1_u2": RunConfig(
            name="mixed_u1_u2",
            lattice_type="kagome",
            rt=1.0,
            ru1=-0.4,
            ru2=1.6,
            nbos=6,
            nlx=2,
            nly=2,
            ltrot=10,
            beta=2.0,
            nwrap=2,
            nbin=400,
            nsweep=8,
            nthermal=200,
        ),
        "strong_u2": RunConfig(
            name="strong_u2",
            lattice_type="kagome",
            rt=1.0,
            ru1=0.0,
            ru2=3.0,
            nbos=6,
            nlx=2,
            nly=2,
            ltrot=12,
            beta=2.4,
            nwrap=3,
            nbin=400,
            nsweep=8,
            nthermal=200,
        ),
        "mixed_nbos3": RunConfig(
            name="mixed_nbos3",
            lattice_type="kagome",
            rt=1.0,
            ru1=-0.4,
            ru2=1.6,
            nbos=3,
            nlx=2,
            nly=2,
            ltrot=10,
            beta=2.0,
            nwrap=2,
            nbin=400,
            nsweep=8,
            nthermal=200,
        ),
        "mixed_l3x2_nbos9": RunConfig(
            name="mixed_l3x2_nbos9",
            lattice_type="kagome",
            rt=1.0,
            ru1=-0.4,
            ru2=1.6,
            nbos=9,
            nlx=3,
            nly=2,
            ltrot=10,
            beta=2.0,
            nwrap=2,
            nbin=2400,
            nsweep=8,
            nthermal=1200,
        ),
        "triangular_weak_u2": RunConfig(
            name="triangular_weak_u2",
            lattice_type="triangular",
            rt=1.0,
            ru1=0.0,
            ru2=1.0,
            nbos=6,
            nlx=2,
            nly=2,
            ltrot=8,
            beta=1.6,
            nwrap=2,
            nbin=120,
            nsweep=8,
            nthermal=60,
            is_warm=True,
            nwarm=20,
            ini_ham=5,
            ini_twist=1.0e-4,
        ),
    }


def lattice_norb(lattice_type: str) -> int:
    key = lattice_type.strip().lower()
    if key == "triangular":
        return 1
    if key == "kagome":
        return 3
    raise ValueError(f"Unsupported lattice_type for warm-start validation: {lattice_type}")


def expected_confin_lines(cfg: RunConfig) -> int:
    ndim = cfg.nlx * cfg.nly * lattice_norb(cfg.lattice_type)
    return 1 + N_AUX * ndim * cfg.ltrot


def validate_confin_file(path: Path, cfg: RunConfig) -> tuple[bool, str]:
    if not path.exists():
        return False, "missing file"
    with path.open("r", encoding="ascii") as stream:
        line_count = sum(1 for _ in stream)
    expected = expected_confin_lines(cfg)
    if line_count != expected:
        return False, f"expected {expected} lines, found {line_count}"
    return True, ""


def fortran_bool(value: bool) -> str:
    return ".true." if value else ".false."


def production_parameter_sets(
    *,
    lattice_type: str = "triangular",
    l_values: tuple[int, ...] = (12,),
    nbos_values: tuple[int, ...] = (100000, 1000000, 10000000),
    u2_values: tuple[float, ...] = (100.0, 1000.0, 10000.0),
    rt: float = 1.0,
    ru1: float = 0.0,
    beta: float = 256.0,
    dtau: float = 1.0e-3,
    nwrap: int = 32,
    nbin: int = 64,
    nsweep: int = 1,
    nthermal: int = 32,
    ini_type: int = 2,
    ini_ampl: float = 0.1,
    ini_type_values: tuple[int, ...] | None = None,
    ini_ampl_values: tuple[float, ...] | None = None,
    ini_ham: int = 5,
    ini_twist: float = 1.0e-4,
    imbalance: float = 0.0,
) -> dict[str, RunConfig]:
    ltrot = int(round(beta / dtau))
    configs: dict[str, RunConfig] = {}
    ini_types = ini_type_values if ini_type_values else (ini_type,)
    ini_ampls = ini_ampl_values if ini_ampl_values else (ini_ampl,)
    multi_init = len(ini_types) > 1 or len(ini_ampls) > 1

    def format_ampl_tag(value: float) -> str:
        return f"{value:.6g}".replace("+", "").replace("-", "m").replace(".", "p")

    for lval in l_values:
        for nbos in nbos_values:
            for u2 in u2_values:
                base_name = f"{lattice_type}_L{lval}_N{nbos}_U2_{u2:.0f}".replace(".", "p")
                for ini_type_item in ini_types:
                    for ini_ampl_item in ini_ampls:
                        name = base_name
                        if multi_init:
                            name = f"{base_name}_iniT{ini_type_item}_A{format_ampl_tag(ini_ampl_item)}"
                        configs[name] = RunConfig(
                            name=name,
                            lattice_type=lattice_type,
                            rt=rt,
                            ru1=ru1,
                            ru2=u2,
                            nbos=nbos,
                            nlx=lval,
                            nly=lval,
                            ltrot=ltrot,
                            beta=beta,
                            nwrap=nwrap,
                            nbin=nbin,
                            nsweep=nsweep,
                            nthermal=nthermal,
                            ini_type=ini_type_item,
                            ini_ampl=ini_ampl_item,
                            ini_ham=ini_ham,
                            ini_twist=ini_twist,
                            imbalance=imbalance,
                        )
    return configs


def write_param_file(
    path: Path,
    cfg: RunConfig,
    *,
    is_global: bool,
    nfrog: int,
    hmc_dt: float,
    hmc_jitter: int = 0,
    hmc_mass: float = 1.0,
    hmc_block_tau: int = 0,
    hmc_block_sites: int = 0,
) -> None:
    nlx_therm, nly_therm, ltrot_therm = cfg.therm_dims()
    lines = [
        cfg.lattice_type,
        f"{cfg.rt} {cfg.ru1} {cfg.ru2} {cfg.nbos}",
        f"{cfg.nlx} {cfg.nly} {cfg.ltrot} {cfg.beta}",
        f"{nlx_therm} {nly_therm} {ltrot_therm}",
        f"{cfg.nwrap} {cfg.nbin} {cfg.nsweep} {cfg.shift_loc}",
        f"{fortran_bool(cfg.is_tau)} {cfg.nthermal}",
        f"{fortran_bool(cfg.is_warm)} {cfg.nwarm} {cfg.shift_warm_1} {cfg.shift_warm_2}",
        f"{fortran_bool(is_global)} {nfrog} {hmc_dt} {hmc_jitter} {hmc_mass} {hmc_block_tau} {hmc_block_sites}",
        f"{cfg.ini_type} {cfg.ini_ampl} {cfg.ini_bias_1} {cfg.ini_bias_2}",
        f"{cfg.ini_ham} {cfg.ini_twist} {cfg.imbalance}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_seed_files(run_dir: Path, seed: int, *, confin_from: Path | None = None) -> None:
    if confin_from is None:
        (run_dir / "confin.txt").write_text("0\n", encoding="ascii")
    else:
        lines = confin_from.read_text(encoding="ascii").splitlines()
        if not lines:
            raise ValueError(f"Warm-start file is empty: {confin_from}")
        lines[0] = str(seed)
        (run_dir / "confin.txt").write_text("\n".join(lines) + "\n", encoding="ascii")
    seeds = [seed + offset for offset in range(8)]
    (run_dir / "seeds.txt").write_text("\n".join(str(item) for item in seeds) + "\n", encoding="ascii")


def prepare_run_dir(
    run_dir: Path,
    cfg: RunConfig,
    *,
    is_global: bool,
    nfrog: int,
    hmc_dt: float,
    seed: int,
    hmc_jitter: int = 0,
    hmc_mass: float = 1.0,
    hmc_block_tau: int = 0,
    hmc_block_sites: int = 0,
    binary: Path = DEFAULT_BINARY,
    confin_from: Path | None = None,
) -> None:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    shutil.copy2(binary, run_dir / "BPQMC.out")
    write_seed_files(run_dir, seed, confin_from=confin_from)
    write_param_file(
        run_dir / "paramC_sets.txt",
        cfg,
        is_global=is_global,
        nfrog=nfrog,
        hmc_dt=hmc_dt,
        hmc_jitter=hmc_jitter,
        hmc_mass=hmc_mass,
        hmc_block_tau=hmc_block_tau,
        hmc_block_sites=hmc_block_sites,
    )


def run_case(run_dir: Path, *, np_ranks: int = 1, env_overrides: dict[str, str] | None = None) -> None:
    env = os.environ.copy()
    if env_overrides:
        env.update(env_overrides)
    with (run_dir / "output.log").open("w", encoding="ascii") as stream:
        subprocess.run(
            ["mpirun", "-np", str(np_ranks), "./BPQMC.out"],
            cwd=run_dir,
            stdout=stream,
            stderr=subprocess.STDOUT,
            env=env,
            check=True,
        )


def parse_info_metrics(run_dir: Path) -> dict[str, float]:
    metrics: dict[str, float] = {}
    info_path = run_dir / "info.txt"
    for line in info_path.read_text(encoding="ascii").splitlines():
        match = ACCEPT_RE.match(line)
        if match:
            key = "_".join(match.group(1).split())
            metrics[key] = float(match.group(2))
    return metrics


def read_scalar_series(run_dir: Path, name: str) -> list[float]:
    path = run_dir / name
    if not path.exists():
        raise FileNotFoundError(path)
    values = []
    for line in path.read_text(encoding="ascii").splitlines():
        stripped = line.strip()
        if stripped:
            values.append(float(stripped.split()[0]))
    return values


def read_complex_series_real(run_dir: Path, name: str) -> list[float]:
    path = run_dir / name
    if not path.exists():
        raise FileNotFoundError(path)
    values = []
    for line in path.read_text(encoding="ascii").splitlines():
        stripped = line.strip()
        if stripped:
            values.append(float(stripped.split()[0]))
    return values


def series_mean(values: list[float]) -> float:
    return mean(values)


def sample_stderr(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return 0.0
    avg = mean(values)
    var = sum((value - avg) ** 2 for value in values) / (n - 1)
    return math.sqrt(var / n)


def lag1_autocorr(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return 0.0
    avg = mean(values)
    centered = [value - avg for value in values]
    denom = sum(value * value for value in centered)
    if denom == 0.0:
        return 0.0
    numer = sum(centered[idx] * centered[idx + 1] for idx in range(n - 1))
    return numer / denom


def integrated_autocorr_time(values: list[float], max_lag: int | None = None) -> float:
    n = len(values)
    if n < 2:
        return 0.5
    avg = mean(values)
    centered = [value - avg for value in values]
    var = sum(value * value for value in centered) / n
    if var == 0.0:
        return 0.5
    limit = n // 2 if max_lag is None else min(max_lag, n - 1)
    tau = 0.5
    for lag in range(1, limit + 1):
        cov = sum(centered[idx] * centered[idx + lag] for idx in range(n - lag)) / (n - lag)
        rho = cov / var
        if lag > 1 and rho <= 0.0:
            break
        tau += rho
    return max(tau, 0.5)


def effective_sample_size(values: list[float]) -> float:
    if not values:
        return 0.0
    span = max(values) - min(values)
    scale = max(1.0, max(abs(value) for value in values))
    if span <= 1.0e-12 * scale:
        return 0.0
    tau = integrated_autocorr_time(values)
    return len(values) / (2.0 * tau)


def ess_per_second(values: list[float], total_cpu_time: float) -> float:
    if total_cpu_time <= 0.0:
        return 0.0
    return effective_sample_size(values) / total_cpu_time
