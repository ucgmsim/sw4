"""Test case specifications (tests/cases/<group>/<name>.toml)."""

import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

TIERS = ("pr", "nightly")


@dataclass
class Run:
    name: str
    tier: str
    params: dict[str, Any]
    ranks: int
    threads: int
    timeout: int
    expect: str  # "success" | "failure"
    stdout_regex: str | None
    input: str
    generate: str | None
    assets: list[str]
    xfail: str | dict | None = None  # known defect (or {precision: reason}): expectation currently NOT met
    after: list[str] = field(default_factory=list)  # runs whose outputs this run reads
    copy_from: dict = field(default_factory=dict)  # {run = name, files = [...]}: copied in before running


@dataclass
class Check:
    name: str
    kind: str
    tier: str
    runs: list[str]
    options: dict[str, Any]
    # Checks that compare against another build's outputs (float vs double).
    needs_peer: bool = False


@dataclass
class Case:
    id: str  # "<group>/<name>"
    path: Path
    description: str
    labels: list[str]
    requires: list[str]
    precisions: list[str]
    runs: dict[str, Run]
    checks: dict[str, Check]
    raw: dict[str, Any] = field(repr=False, default_factory=dict)

    @property
    def dir(self) -> Path:
        return self.path.parent


def _tier(value: str, where: str) -> str:
    if value not in TIERS:
        raise ValueError(f"{where}: tier must be one of {TIERS}, got {value!r}")
    return value


def load_case(path: Path, cases_dir: Path) -> Case:
    raw = tomllib.loads(path.read_text())
    case_id = path.relative_to(cases_dir).with_suffix("").as_posix()
    base_params = raw.get("params", {})
    runs: dict[str, Run] = {}
    for r in raw.get("run", []):
        name = r["name"]
        if name in runs:
            raise ValueError(f"{case_id}: duplicate run {name}")
        runs[name] = Run(
            name=name,
            tier=_tier(r.get("tier", "pr"), f"{case_id}:{name}"),
            params={**base_params, **r.get("params", {})},
            ranks=int(r.get("ranks", raw.get("ranks", 4))),
            threads=int(r.get("threads", raw.get("threads", 1))),
            timeout=int(r.get("timeout", raw.get("timeout", 300))),
            expect=r.get("expect", raw.get("expect", "success")),
            stdout_regex=r.get("stdout_regex", raw.get("stdout_regex")),
            input=r.get("input", raw.get("input", "")),
            generate=r.get("generate", raw.get("generate")),
            assets=list(r.get("assets", raw.get("assets", []))),
            xfail=r.get("xfail"),
            after=list(r.get("after", [])),
            copy_from=dict(r.get("copy_from", {})),
        )
    checks: dict[str, Check] = {}
    for c in raw.get("check", []):
        name = c.get("name", c["kind"])
        if name in checks:
            raise ValueError(f"{case_id}: duplicate check {name}")
        crun = list(c.get("runs", []))
        for rn in crun:
            if rn not in runs:
                raise ValueError(f"{case_id}:{name}: unknown run {rn}")
        default_tier = c.get("tier")
        if default_tier is None:
            used = crun or list(runs)
            default_tier = "nightly" if any(runs[r].tier == "nightly" for r in used) else "pr"
        else:
            # A pr check may only depend on pr runs, else -L pr would pull nightly runs in.
            if default_tier == "pr":
                bad = [r for r in crun if runs[r].tier != "pr"]
                if bad:
                    raise ValueError(f"{case_id}:{name}: pr check depends on nightly runs {bad}")
        opts = {k: v for k, v in c.items() if k not in ("name", "kind", "tier", "runs")}
        checks[name] = Check(
            name=name,
            kind=c["kind"],
            tier=_tier(default_tier, f"{case_id}:{name}"),
            runs=crun,
            options=opts,
            needs_peer=bool(c.get("peer", False)),
        )
    return Case(
        id=case_id,
        path=path,
        description=raw.get("description", "").strip(),
        labels=list(raw.get("labels", [])),
        requires=list(raw.get("requires", [])),
        precisions=list(raw.get("precisions", ["double", "float"])),
        runs=runs,
        checks=checks,
        raw=raw,
    )


def load_cases(cases_dir: Path) -> list[Case]:
    return [load_case(p, cases_dir) for p in sorted(cases_dir.rglob("*.toml"))]


def find_case(cases_dir: Path, case_id: str) -> Case:
    return load_case(cases_dir / f"{case_id}.toml", cases_dir)


_PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z_][A-Za-z0-9_]*)\s*\}\}")


def render(template: str, params: dict[str, Any]) -> str:
    def sub(m: re.Match) -> str:
        key = m.group(1)
        if key not in params:
            raise KeyError(f"template placeholder {{{{{key}}}}} has no value")
        v = params[key]
        if isinstance(v, bool):
            return "1" if v else "0"
        if isinstance(v, float):
            return repr(v)
        return str(v)

    return _PLACEHOLDER.sub(sub, template)
