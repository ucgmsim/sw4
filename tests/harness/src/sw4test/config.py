"""Build configuration shared by every harness command (written by `configure`)."""

import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class Config:
    sw4: str
    precision: str  # "double" | "float"
    features: list[str]  # subset of {"hdf5", "proj", "fftw"}
    source_dir: str
    out_dir: str
    mpiexec: str
    numproc_flag: str = "-n"
    mpi_preflags: list[str] = field(default_factory=list)
    uv: str = "uv"
    # Directory of another build's tests/out, for float-vs-double cross checks.
    peer_out_dir: str | None = None
    peer_precision: str | None = None

    @property
    def cases_dir(self) -> Path:
        return Path(self.source_dir) / "tests" / "cases"

    @property
    def golden_dir(self) -> Path:
        return Path(self.source_dir) / "tests" / "golden" / self.precision

    def golden_dir_for(self, precision: str) -> Path:
        return Path(self.source_dir) / "tests" / "golden" / precision

    def run_dir(self, case: str, run: str, out_dir: str | None = None) -> Path:
        return Path(out_dir or self.out_dir) / case / run

    def mpi_command(self, ranks: int) -> list[str]:
        return [self.mpiexec, self.numproc_flag, str(ranks), *self.mpi_preflags]

    def save(self, path: Path) -> None:
        path.write_text(json.dumps(asdict(self), indent=2) + "\n")

    @classmethod
    def load(cls, path: str | Path) -> "Config":
        return cls(**json.loads(Path(path).read_text()))


def detect_mpi_preflags(mpiexec: str) -> list[str]:
    """Flags so several MPI tests can share a machine without fighting.

    Open MPI binds each rank to a core by default, counting from core 0, so
    concurrent ctest jobs would all pile onto the same cores. It also refuses
    more ranks than it thinks there are slots (a 4-vCPU CI runner reports 2).
    """
    exe = shutil.which(mpiexec) or mpiexec
    try:
        out = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.TimeoutExpired):
        return []
    if "Open MPI" in out or "OpenRTE" in out:
        return ["--oversubscribe", "--bind-to", "none"]
    if re.search(r"HYDRA|MPICH", out):
        return ["-bind-to", "none"]
    return []
