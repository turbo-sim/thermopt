"""Keep project-owned property access behind jaxprop, including study scripts."""

import ast
import json
import os
from pathlib import Path
import tomllib


ROOT = Path(__file__).resolve().parents[1]


def project_sources():
    excluded = {".git", ".venv", "venv", "__pycache__", "node_modules", ".pytest_cache", ".test_artifacts", "_build"}
    for directory, dirs, files in os.walk(ROOT):
        dirs[:] = [name for name in dirs if name not in excluded and not name.startswith(".test_tmp")]
        for name in files:
            path = Path(directory) / name
            if path.suffix in {".py", ".ipynb"}:
                yield path


def test_project_property_access_uses_jaxprop():
    violations = []
    for path in project_sources():
        text = path.read_text(encoding="utf-8-sig")
        if path.suffix == ".ipynb":
            text = "\n".join("".join(cell["source"]) for cell in json.loads(text)["cells"] if cell["cell_type"] == "code")
        for node in ast.walk(ast.parse(text, filename=str(path))):
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                modules = [node.module or ""]
            if any(module.split(".")[0].lower() in {"coolprop", "coolpropx"} for module in modules):
                violations.append(f"{path.relative_to(ROOT)}:{node.lineno}: direct backend import")
            if isinstance(node, ast.Attribute) and node.attr in {"abstract_state", "PropsSI", "AbstractState"}:
                violations.append(f"{path.relative_to(ROOT)}:{node.lineno}: raw backend access")
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr in {"import_module", "__import__"}:
                if node.args and isinstance(node.args[0], ast.Constant) and str(node.args[0].value).lower().startswith(("coolprop", "coolpropx")):
                    violations.append(f"{path.relative_to(ROOT)}:{node.lineno}: dynamic backend import")
    assert not violations, "\n".join(violations)


def test_direct_property_dependency_is_jaxprop():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    poetry = config["tool"]["poetry"]
    groups = [poetry["dependencies"], *(group.get("dependencies", {}) for group in poetry.get("group", {}).values())]
    for dependencies in groups:
        assert not {name.lower() for name in dependencies} & {"coolprop", "coolpropx"}
    assert "jaxprop" in poetry["dependencies"]
