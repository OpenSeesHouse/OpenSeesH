#!/usr/bin/env python3
"""Audit Win64 ClCompile vs Make OBJS vs SRC/Makefile *_LIBS.

Reports set differences used before Linux/HPC rebuilds.
Exit 0 if P0/P1 focus diffs are empty (see FOCUS_*); exit 1 otherwise.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "SRC"
WIN64 = ROOT / "Win64" / "proj"
MAKEFILE = SRC / "Makefile"

ALLOW_PROJ_GLOBS = [
    "element/*.vcxproj",
    "material/*.vcxproj",
    "domain/*.vcxproj",
    "recorder/*.vcxproj",
    "openSeesPy/*.vcxproj",
    "modelbuilder/*.vcxproj",
]

# Focus dirs where mismatch blocks Python-only archive/link
FOCUS_PREFIXES = (
    "damping/",
    "material/uniaxial/PY/",
    "element/UP-ucsd/",
    "api/",
    "damage/",
    "recorder/",
    "renderer/",
)

SKIP_NAME_RE = re.compile(r"(?i)(^|/)(Tcl|tclMain|elementAPI_TCL)")


def norm_rel(p: str) -> str:
    p = p.replace("\\", "/")
    # Drop leading ../.. noise from vcxproj
    if "SRC/" in p:
        p = p[p.index("SRC/") + 4 :]
    p = p.replace("Src/", "").replace("src/", "")
    # Case-fix PY materials folder
    p = re.sub(r"(?i)material/uniaxial/py/", "material/uniaxial/PY/", p)
    if p.endswith(".cpp"):
        p = p[:-4] + ".o"
    elif p.endswith(".c"):
        p = p[:-2] + ".o"
    elif p.endswith(".f") or p.endswith(".f90"):
        p = re.sub(r"\.(f90|f)$", ".o", p, flags=re.I)
    return p


def win64_allowlist() -> set[str]:
    out: set[str] = set()
    for glob in ALLOW_PROJ_GLOBS:
        for proj in WIN64.glob(glob):
            text = proj.read_text(encoding="utf-8", errors="ignore")
            for m in re.finditer(r'<ClCompile\s+Include="([^"]+)"', text):
                rel = norm_rel(m.group(1))
                if SKIP_NAME_RE.search(rel):
                    continue
                if not rel.endswith(".o"):
                    continue
                out.add(rel)
    return out


def parse_makefile_objs(makefile: Path, varname: str = "OBJS") -> set[str]:
    if not makefile.is_file():
        return set()
    text = makefile.read_text(encoding="utf-8", errors="ignore")
    # Prefer GRAPHICS=NONE branch when present (Ubuntu24 Python profile)
    none = re.search(
        r"(?ms)^ifeq\s*\(\$\(GRAPHICS\),\s*NONE\)\s*\r?\n(?:\r?\n)*OBJS\s*=\s*((?:.*\\\r?\n)*.*?)(?=\r?\n#endif|\r?\nifeq|\r?\nelse|\Z)",
        text,
    )
    if none and varname == "OBJS":
        body = none.group(1)
        return {
            tok
            for tok in re.findall(r"([\w./+-]+\.o)", body)
            if not SKIP_NAME_RE.search(tok)
        }
    objs: set[str] = set()
    for m in re.finditer(
        rf"(?m)^(?:{varname}|BLOCK_OBJS)\s*\+?=\s*((?:.*\\\n)*.*)", text
    ):
        body = m.group(1)
        for tok in re.findall(r"([\w./+-]+\.o)", body):
            if SKIP_NAME_RE.search(tok):
                continue
            objs.add(tok)
    return objs


def walk_active_built() -> set[str]:
    """Approximate active Make tree for Ubuntu24 PYTHON (no reliability/tcl)."""
    built: set[str] = set()

    def add_dir(rel_dir: str, objs: set[str]) -> None:
        for o in objs:
            built.add(f"{rel_dir}/{o}".replace("//", "/"))

    top = [
        "matrix",
        "database",
        "system_of_eqn",
        "modelbuilder",
        "tagged",
        "material",
        "analysis",
        "graph",
        "utility",
        "convergenceTest",
        "recorder",
        "renderer",
        "domain",
        "handler",
        "actor",
        "machine",
        "element",
        "damping",
        "coordTransformation",
        "api",
        "damage",
        "string",
    ]
    for d in top:
        mf = SRC / d / "Makefile"
        objs = parse_makefile_objs(mf)
        if d == "modelbuilder":
            add_dir(d, objs)
            text = (SRC / "modelbuilder" / "tcl" / "Makefile").read_text(
                encoding="utf-8", errors="ignore"
            )
            bm = re.search(r"(?m)^BLOCK_OBJS\s*=\s*((?:.*\\\n)*.*)", text)
            if bm:
                for tok in re.findall(r"([\w./+-]+\.o)", bm.group(1)):
                    built.add(f"modelbuilder/tcl/{tok}")
            continue
        add_dir(d, objs)

    ele_mf = (SRC / "element" / "Makefile").read_text(encoding="utf-8", errors="ignore")
    for m in re.finditer(r"\$\(FE\)/element/([\w+-]+)", ele_mf):
        sub = m.group(1)
        objs = parse_makefile_objs(SRC / "element" / sub / "Makefile")
        add_dir(f"element/{sub}", objs)
        sub_mf = SRC / "element" / sub / "Makefile"
        if sub_mf.is_file():
            t = sub_mf.read_text(encoding="utf-8", errors="ignore")
            for n in re.finditer(
                r"\$\(CD\)\s+\$\(FE\)/element/" + re.escape(sub) + r"/([\w+-]+)", t
            ):
                nest = n.group(1)
                nobjs = parse_makefile_objs(SRC / "element" / sub / nest / "Makefile")
                add_dir(f"element/{sub}/{nest}", nobjs)

    for sub in ("uniaxial", "nD", "section"):
        objs = parse_makefile_objs(SRC / "material" / sub / "Makefile")
        add_dir(f"material/{sub}", objs)
        mf = SRC / "material" / sub / "Makefile"
        if not mf.is_file():
            continue
        t = mf.read_text(encoding="utf-8", errors="ignore")
        for n in re.finditer(rf"\$\(CD\)\s+\$\(FE\)/material/{sub}/([\w+-]+)", t):
            nest = n.group(1)
            nobjs = parse_makefile_objs(SRC / "material" / sub / nest / "Makefile")
            add_dir(f"material/{sub}/{nest}", nobjs)
            nest_mf = SRC / "material" / sub / nest / "Makefile"
            if nest_mf.is_file():
                nt = nest_mf.read_text(encoding="utf-8", errors="ignore")
                for n2 in re.finditer(
                    rf"\$\(CD\)\s+\$\(FE\)/material/{sub}/{nest}/([\w+-]+)", nt
                ):
                    nest2 = n2.group(1)
                    n2objs = parse_makefile_objs(
                        SRC / "material" / sub / nest / nest2 / "Makefile"
                    )
                    add_dir(f"material/{sub}/{nest}/{nest2}", n2objs)

    robjs = parse_makefile_objs(SRC / "recorder" / "response" / "Makefile")
    add_dir("recorder/response", robjs)
    return built


def archived_libs(graphics: str = "NONE") -> set[str]:
    text = MAKEFILE.read_text(encoding="utf-8", errors="ignore")
    text_no_tcl = re.sub(
        r"(?ms)^TCL_COMMAND_LIBS\s*=.*?(?=^[A-Z][A-Z0-9_]*\s*=)",
        "",
        text,
    )
    # Keep only GRAPHICS=NONE RENDERER_LIBS assignment
    if graphics == "NONE":
        text_no_tcl = re.sub(
            r"(?ms)^ifeq\s*\(\$\(GRAPHICS\),\s*UsingOpenGL\).*?^endif\s*\n",
            "",
            text_no_tcl,
        )
        # Remove the else X11 RENDERER block already consumed; keep NONE block
        text_no_tcl = re.sub(
            r"(?ms)^RENDERER_LIBS\s*=(?!.*PlainMap)(?:.*\\\n)*.*?\n\n",
            "",
            text_no_tcl,
        )
    archived: set[str] = set()
    for m in re.finditer(r"\$\(FE\)/([^\s\\]+\.o)", text_no_tcl):
        rel = norm_rel(m.group(1))
        if SKIP_NAME_RE.search(rel):
            continue
        # FilePlotter only when graphics enabled
        if graphics == "NONE" and rel.endswith("recorder/FilePlotter.o"):
            continue
        archived.add(rel)
    return archived


def focus(s: set[str]) -> set[str]:
    return {x for x in s if x.startswith(FOCUS_PREFIXES)}


def main() -> int:
    allow = win64_allowlist()
    built = walk_active_built()
    # OpsRecorderCommands is built via interpreter Makefile, not recorder OBJS
    built.add("recorder/OpsRecorderCommands.o")
    archived = archived_libs("NONE")

    built_m_allow = focus(built - allow)
    allow_m_built = focus(allow - built)
    arch_m_built = focus(archived - built)
    built_m_arch = focus(built - archived)

    print("=== FOCUS audit (damping/PY/UP-ucsd/api/damage/recorder/renderer) ===")
    print(f"allowlist total: {len(allow)}  built: {len(built)}  archived: {len(archived)}")
    for title, s in [
        ("Built - Allowlist (extra in Make)", built_m_allow),
        ("Allowlist - Built (missing from Make)", allow_m_built),
        ("Archived - Built (ar will fail / stale)", arch_m_built),
        ("Built - Archived (silent drop from lib)", built_m_arch),
    ]:
        print(f"\n-- {title}: {len(s)}")
        for x in sorted(s)[:80]:
            print(" ", x)
        if len(s) > 80:
            print(f"  ... +{len(s) - 80} more")

    hard = arch_m_built | {
        x
        for x in built_m_arch
        if x != "recorder/OpsRecorderCommands.o"
        and not x.startswith("renderer/")
        or x in ("renderer/PlainMap.o", "renderer/Renderer.o")
        and x in built_m_arch
    }
    # OpsRecorderCommands links via interpreter, not FE_LIBRARY.
    # Renderer test mains (main.o, …) are not production OBJS under GRAPHICS=NONE.
    hard = set(arch_m_built)
    for x in built_m_arch:
        if x == "recorder/OpsRecorderCommands.o":
            continue
        if x.startswith("renderer/") and x not in (
            "renderer/PlainMap.o",
            "renderer/Renderer.o",
        ):
            continue
        hard.add(x)
    if hard:
        print("\nAUDIT_FAIL focus archive/build mismatches")
        return 1
    print("\nAUDIT_OK focus Built<->Archived")
    return 0


if __name__ == "__main__":
    sys.exit(main())
