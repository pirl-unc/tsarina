"""Exercise development installation with real pip and an offline local corpus."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import venv
import zipfile
from pathlib import Path

import pytest

REPOSITORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("develop", REPOSITORY / "scripts/develop.py")
develop = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(develop)

# A stdlib-only PEP 660 backend keeps the resolver tests offline and independent
# of whichever setuptools/wheel versions happen to be installed in CI.
BACKEND = """
import json
import pathlib
import zipfile

def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    root = pathlib.Path(__file__).parent
    data = json.loads((root / "project.json").read_text())
    name, version = data["name"], data["version"]
    dist_info = f"{name}-{version}.dist-info"
    wheel_name = f"{name}-{version}-py3-none-any.whl"
    fields = ["Metadata-Version: 2.1", f"Name: {name}", f"Version: {version}",
              "Provides-Extra: dev"]
    fields += [f"Requires-Dist: {r}" for r in data["requirements"]]
    with zipfile.ZipFile(pathlib.Path(wheel_directory) / wheel_name, "w") as wheel:
        wheel.writestr(f"{dist_info}/METADATA", "\\n".join(fields) + "\\n")
        wheel.writestr(f"{dist_info}/WHEEL",
                      "Wheel-Version: 1.0\\nRoot-Is-Purelib: true\\nTag: py3-none-any\\n")
        wheel.writestr(f"{name}.pth", str(root) + "\\n")
        wheel.writestr(f"{dist_info}/RECORD", "")
    return wheel_name
"""


def _wheel(directory, name, version, requirements=()):
    path = directory / f"{name}-{version}-py3-none-any.whl"
    dist_info = f"{name}-{version}.dist-info"
    fields = ["Metadata-Version: 2.1", f"Name: {name}", f"Version: {version}"]
    fields += [f"Requires-Dist: {requirement}" for requirement in requirements]
    with zipfile.ZipFile(path, "w") as wheel:
        wheel.writestr(f"{dist_info}/METADATA", "\n".join(fields) + "\n")
        wheel.writestr(
            f"{dist_info}/WHEEL",
            "Wheel-Version: 1.0\nRoot-Is-Purelib: true\nTag: py3-none-any\n",
        )
        wheel.writestr(f"{name}/__init__.py", f"__version__ = {version!r}\n")
        wheel.writestr(f"{dist_info}/RECORD", "")
    return path


def _project(directory, name, version, requirements=()):
    root = directory / name
    root.mkdir()
    (root / "pyproject.toml").write_text(
        '[build-system]\nrequires = []\nbuild-backend = "backend"\nbackend-path = ["."]\n'
    )
    (root / "backend.py").write_text(BACKEND)
    (root / "project.json").write_text(
        json.dumps({"name": name, "version": version, "requirements": list(requirements)})
    )
    (root / name).mkdir()
    (root / name / "__init__.py").write_text(f"__version__ = {version!r}\n")
    return root


@pytest.fixture
def installation(tmp_path):
    environment = tmp_path / "venv"
    venv.EnvBuilder(with_pip=True, symlinks=True).create(environment)
    python = environment / "bin/python"
    wheels = tmp_path / "wheels"
    wheels.mkdir()
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.update(
        VIRTUAL_ENV=str(environment),
        PATH=f"{environment / 'bin'}:/usr/bin:/bin",
        PIP_NO_INDEX="1",
        PIP_FIND_LINKS=str(wheels),
        PIP_CONFIG_FILE=os.devnull,
        PIP_DISABLE_PIP_VERSION_CHECK="1",
        PIP_CACHE_DIR=str(tmp_path / "pip-cache"),
        SIBLING_ROOT=str(tmp_path),
    )
    root = _project(tmp_path, "tsarina", "1.0", ["hitlist>=1"])
    shutil.copy(REPOSITORY / "develop.sh", root)
    (root / "scripts").mkdir()
    shutil.copy(REPOSITORY / "scripts/develop.py", root / "scripts")

    def run(*arguments, check=False):
        return subprocess.run(
            [str(python), *arguments],
            cwd=tmp_path,
            env=env,
            capture_output=True,
            text=True,
            check=check,
        )

    def install(*paths):
        run("-m", "pip", "install", *(str(path) for path in paths), check=True)

    def snapshot():
        return json.loads(
            run(
                "-c",
                "import importlib.metadata as m, json; "
                "print(json.dumps({d.metadata['Name']: d.version for d in m.distributions()}))",
                check=True,
            ).stdout
        )

    def develop_install(*, outside=False):
        return subprocess.run(
            ["/bin/bash", str(root / "develop.sh")],
            cwd=tmp_path if outside else root,
            env=env,
            capture_output=True,
            text=True,
        )

    return wheels, root, run, install, snapshot, develop_install


def test_stale_sibling_is_rejected_before_changing_installed_packages(installation, tmp_path):
    wheels, _, run, install, snapshot, execute = installation
    install(_wheel(wheels, "mhcgnomes", "3.64.4"))
    _wheel(wheels, "hitlist", "1.0", ["mhcgnomes>=3.64.4"])
    _project(tmp_path, "hitlist", "1.0", ["mhcgnomes>=3.64.4"])
    _project(tmp_path, "mhcgnomes", "3.64.2")
    before = snapshot()

    result = execute()

    assert result.returncode != 0, result.stdout + result.stderr
    assert "mhcgnomes" in result.stderr
    assert "3.64.4" in result.stdout + result.stderr
    assert "installed packages were not changed" in result.stderr
    assert snapshot() == before
    assert run("-m", "pip", "check").returncode == 0


def test_stale_sibling_cannot_break_an_untouched_installed_consumer(installation, tmp_path):
    wheels, _, run, install, snapshot, execute = installation
    install(
        _wheel(wheels, "mhcgnomes", "3.64.4"),
        _wheel(wheels, "vaxrank", "3.33.0", ["mhcgnomes>=3.64.4"]),
    )
    _wheel(wheels, "hitlist", "1.0", ["mhcgnomes>=3.20"])
    _project(tmp_path, "hitlist", "1.0", ["mhcgnomes>=3.20"])
    _project(tmp_path, "mhcgnomes", "3.64.2")
    before = snapshot()

    result = execute()

    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "vaxrank 3.33.0 has requirement mhcgnomes>=3.64.4" in output
    assert "3.64.2" in output
    assert "Update the incompatible sibling checkout" in output
    assert "separate virtualenv" in output
    assert snapshot() == before
    assert run("-m", "pip", "check").returncode == 0


def test_joint_install_adds_new_dependency_and_preserves_editables(installation, tmp_path):
    wheels, root, run, install, _, execute = installation
    install(_wheel(wheels, "hitlist", "1.0"), _wheel(wheels, "mhcgnomes", "3.64.2"))
    _wheel(wheels, "new_dependency", "2.0")
    hitlist = _project(tmp_path, "hitlist", "1.1", ["mhcgnomes>=3.64.4", "new-dependency>=2"])
    mhcgnomes = _project(tmp_path, "mhcgnomes", "3.64.4")

    result = execute()

    assert result.returncode == 0, result.stdout + result.stderr
    check = run("-m", "pip", "check")
    assert check.returncode == 0, check.stdout + check.stderr
    for name, path in [("tsarina", root), ("hitlist", hitlist), ("mhcgnomes", mhcgnomes)]:
        installed = json.loads(
            run(
                "-c",
                f"from importlib.metadata import distribution; "
                f"print(distribution({name!r}).read_text('direct_url.json'))",
                check=True,
            ).stdout
        )
        assert installed["dir_info"]["editable"]
        assert installed["url"] == path.as_uri()
        import_path = run("-c", f"import {name}; print({name}.__file__)", check=True).stdout.strip()
        assert Path(import_path) == path / name / "__init__.py"
    assert (
        run("-c", "import new_dependency; print(new_dependency.__version__)").stdout.strip()
        == "2.0"
    )


def test_absent_sibling_uses_release_dependency(installation, tmp_path):
    wheels, _, run, _, _, execute = installation
    _wheel(wheels, "hitlist", "1.0", ["new-dependency>=2; python_version >= '3'"])
    _wheel(wheels, "new_dependency", "2.0")
    _wheel(wheels, "irrelevant_dependency", "1.0")
    _project(tmp_path, "mhcgnomes", "3.64.4", ["irrelevant-dependency; python_version < '2'"])

    result = execute(outside=True)

    assert result.returncode == 0, result.stdout + result.stderr
    assert run("-c", "import hitlist; print(hitlist.__version__)").stdout.strip() == "1.0"
    assert run("-m", "pip", "check").returncode == 0
    assert run("-c", "import irrelevant_dependency").returncode != 0


def test_final_dependency_check_failure_is_reported(monkeypatch, tmp_path, capsys):
    report = {"version": "1", "install": []}
    monkeypatch.setattr(develop, "_check_projection", lambda *_: None)

    def run(command, **_):
        if "--report" in command:
            Path(command[command.index("--report") + 1]).write_text(json.dumps(report))
        if command[-1] == "check":
            raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(develop.subprocess, "run", run)

    assert develop.main(["-e", str(tmp_path)]) == 1
    assert "Development install failed" in capsys.readouterr().err


def test_unsupported_report_fails_closed():
    with pytest.raises(ValueError, match="Unsupported pip installation report version"):
        develop._projected_metadata({"version": "2", "install": []})


def test_projection_cannot_be_shadowed_by_pythonpath(monkeypatch, tmp_path):
    foreign = tmp_path / "foreign"
    dist_info = foreign / "mhcgnomes-3.64.4.dist-info"
    dist_info.mkdir(parents=True)
    (dist_info / "METADATA").write_text("Metadata-Version: 2.1\nName: mhcgnomes\nVersion: 3.64.4\n")
    monkeypatch.setenv("PYTHONPATH", str(foreign))
    proposed = {
        "mhcgnomes": develop._package_metadata("mhcgnomes", "3.64.2", []),
        "vaxrank": develop._package_metadata("vaxrank", "3.33.0", ["mhcgnomes>=3.64.4"]),
    }

    with pytest.raises(subprocess.CalledProcessError):
        develop._check_projection(proposed, tmp_path)
