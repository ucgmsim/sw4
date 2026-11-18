"""Narrow detected compiler externals to the languages SW4 actually needs.

Run as ``spack python narrow-compiler-languages.py <packages.yaml> [...]``.

Why this exists: ``spack compiler find`` records every language it can detect.
On a host that has ``gdc`` installed it writes

    - spec: gcc@16.1.1 languages:='c,c++,d,fortran'

and Spack then cannot use that external at all -- concretizing against
``%c=gcc@16.1.1`` fails with a bare "Cannot satisfy 'gcc@16.1.1'". The same entry
with ``languages:='c,c++,fortran'`` concretizes fine. SW4 only needs C, C++ and
Fortran, so drop everything else from both the spec's ``languages`` list and the
``extra_attributes:compilers`` mapping.

Rewrites the given files in place; a file with nothing to change is left alone.
"""

import re
import sys

import spack.util.spack_yaml as syaml

# Language names as they appear in the `languages` variant, and the corresponding
# keys under extra_attributes:compilers.
KEEP_LANGUAGES = ("c", "c++", "fortran")
KEEP_COMPILER_KEYS = ("c", "cxx", "fortran")

LANGUAGES_RE = re.compile(r"languages:=('|\")?(?P<langs>[^'\"\s]+)('|\")?")


def narrow_spec(spec: str):
    """Return (new_spec, kept_languages) or (spec, None) if nothing to do."""
    match = LANGUAGES_RE.search(spec)
    if not match:
        return spec, None

    langs = [lang for lang in match.group("langs").split(",") if lang]
    kept = [lang for lang in langs if lang in KEEP_LANGUAGES]
    if kept == langs:
        return spec, None
    if not kept:
        # Not a compiler we can use; leave it for Spack to reject with its own
        # error rather than silently rewriting it into something meaningless.
        return spec, None

    replacement = "languages:='{0}'".format(",".join(kept))
    return spec[: match.start()] + replacement + spec[match.end() :], kept


def narrow_external(external) -> bool:
    spec = external.get("spec", "")
    new_spec, kept = narrow_spec(spec)
    if kept is None:
        return False

    external["spec"] = new_spec

    compilers = external.get("extra_attributes", {}).get("compilers")
    if compilers:
        for key in list(compilers):
            if key not in KEEP_COMPILER_KEYS:
                del compilers[key]

    return True


def narrow_file(path: str) -> int:
    try:
        with open(path) as handle:
            data = syaml.load_config(handle)
    except (OSError, syaml.SpackYAMLError):
        return 0

    if not data:
        return 0

    changed = 0
    for config in (data.get("packages") or {}).values():
        if not isinstance(config, dict):
            continue
        for external in config.get("externals") or []:
            if narrow_external(external):
                changed += 1

    if changed:
        with open(path, "w") as handle:
            syaml.dump_config(data, handle, default_flow_style=False)

    return changed


def main(argv) -> int:
    if not argv:
        print(__doc__.splitlines()[2], file=sys.stderr)
        return 2

    total = sum(narrow_file(path) for path in argv)
    print(total)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
