#!/usr/bin/env python3
import argparse
import ast
import json
import sys
from pathlib import Path


def _strip_magics(code: str) -> str:
    lines = []
    for line in code.splitlines():
        stripped = line.lstrip()
        if stripped.startswith(("%%", "%", "!")):
            continue
        lines.append(line)
    return "\n".join(lines)


def load_code(path: Path) -> str:
    if path.suffix == ".ipynb":
        nb = json.loads(path.read_text())
        cells = []
        for cell in nb.get("cells", []):
            if cell.get("cell_type") != "code":
                continue
            src = "".join(cell.get("source", []))
            cells.append(_strip_magics(src))
        return "\n".join(cells)
    return path.read_text()


def stdlib_modules() -> set:
    names = set()
    try:
        names.update(sys.stdlib_module_names)
    except AttributeError:
        pass
    names.update(
        {
            "__future__",
            "typing",
            "pathlib",
            "sys",
            "os",
            "math",
            "json",
            "re",
            "time",
            "types",
            "fnmatch",
            "itertools",
            "collections",
            "subprocess",
            "logging",
            "argparse",
        }
    )
    return names


def extract_modules(code: str) -> list:
    tree = ast.parse(code)
    modules = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".")[0]
                modules.add(root)
        elif isinstance(node, ast.ImportFrom):
            if node.level and node.level > 0:
                continue
            if node.module:
                root = node.module.split(".")[0]
                modules.add(root)
    return sorted(modules)


def build_script(packages: list) -> str:
    pkgs = "\n  ".join(packages)
    return f"""#!/usr/bin/env bash
set -euo pipefail

packages=(
  {pkgs}
)

missing=()
while IFS= read -r pkg; do
  [ -z "${{pkg}}" ] && continue
  missing+=("${{pkg}}")
done < <(python3 - <<'PY'
import importlib.util

packages = [
{",\n".join([f"    {p!r}" for p in packages])}
]

for pkg in packages:
    if importlib.util.find_spec(pkg) is None:
        print(pkg)
PY
)

if [ ${{#missing[@]}} -eq 0 ]; then
  echo "All packages already installed."
  exit 0
fi

echo "Installing missing packages: ${{missing[*]}}"
python3 -m pip install --user "${{missing[@]}}"
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a pip install-check script from a .py or .ipynb file."
    )
    parser.add_argument("input", type=Path, help="Path to .py or .ipynb file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help='Output script path (default: "install_deps_<input_stem>.sh")',
    )
    args = parser.parse_args()

    code = load_code(args.input)
    modules = extract_modules(code)
    stdlib = stdlib_modules()
    packages = [m for m in modules if m not in stdlib]
    script = build_script(packages)

    output = args.output
    if output is None:
        output = Path(f"install_deps_{args.input.stem}.sh")
    output.write_text(script)
    print(f"done. please run: bash {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
