"""Execute one SW4 run of a case."""

import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

from . import generators
from .cases import Case, Run, render
from .config import Config


def prepare(cfg: Config, case: Case, run: Run, run_dir: Path) -> Path:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    params = dict(run.params)
    params.setdefault("precision", cfg.precision)
    for asset in run.assets:
        shutil.copy(case.dir / asset, run_dir / Path(asset).name)
    if run.copy_from:
        src = cfg.run_dir(case.id, run.copy_from["run"])
        for name in run.copy_from.get("files", []):
            shutil.copy(src / name, run_dir / name)
    if run.generate:
        extra = generators.generate(run.generate, run_dir, params, case.dir)
        params.update(extra or {})
    if not run.input:
        raise ValueError(f"{case.id}:{run.name}: no input template")
    # String params may themselves contain {{placeholders}} (often generator
    # outputs such as {{src_lat}}); expand them before rendering the input.
    for _ in range(2):
        params = {k: (render(v, params) if isinstance(v, str) and "{{" in v else v) for k, v in params.items()}
    text = render((case.dir / run.input).read_text(), params)
    inp = run_dir / "sw4.in"
    inp.write_text(text)
    (run_dir / "params.json").write_text(json.dumps(params, indent=2, default=str) + "\n")
    return inp


def execute(cfg: Config, case: Case, run: Run) -> int:
    run_dir = cfg.run_dir(case.id, run.name)
    inp = prepare(cfg, case, run, run_dir)
    cmd = [*cfg.mpi_command(run.ranks), cfg.sw4, inp.name]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(run.threads)
    print("cwd:", run_dir)
    print("cmd:", " ".join(cmd), flush=True)
    t0 = time.monotonic()
    with open(run_dir / "sw4.out", "w") as out, open(run_dir / "sw4.err", "w") as err:
        proc = subprocess.run(cmd, cwd=run_dir, stdout=out, stderr=err, env=env)
    wall = time.monotonic() - t0
    stdout = (run_dir / "sw4.out").read_text(errors="replace")
    stderr = (run_dir / "sw4.err").read_text(errors="replace")
    result = {"returncode": proc.returncode, "wall_seconds": round(wall, 3), "ranks": run.ranks,
              "threads": run.threads, "precision": cfg.precision}
    (run_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"sw4 exit code {proc.returncode} after {wall:.1f} s")

    ok = True
    if run.expect == "success":
        if proc.returncode != 0:
            print("FAIL: sw4 exited non-zero; tail of output:")
            print("\n".join((stdout + stderr).splitlines()[-40:]))
            ok = False
    elif run.expect == "failure":
        if proc.returncode == 0:
            print("FAIL: sw4 was expected to reject this input but exited 0")
            ok = False
        if re.search(r"Segmentation|signal 11|SIGSEGV|core dumped", stdout + stderr):
            print("FAIL: expected a clean error, got a crash")
            ok = False
    elif run.expect != "any":
        raise ValueError(f"unknown expect {run.expect!r}")
    if run.stdout_regex and not re.search(run.stdout_regex, stdout + stderr, re.MULTILINE):
        print(f"FAIL: output does not match {run.stdout_regex!r}; tail of output:")
        print("\n".join((stdout + stderr).splitlines()[-40:]))
        ok = False
    xfail = run.xfail.get(cfg.precision) if isinstance(run.xfail, dict) else run.xfail
    if xfail:
        if ok:
            print(f"XPASS: expectations now met; remove xfail from this run ({xfail})")
            return 1
        print(f"XFAIL (known defect, still present): {xfail}")
        return 0
    return 0 if ok else 1
