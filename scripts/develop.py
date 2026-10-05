"""Preflight a joint editable install against the entire selected environment.

Only pip's public CLI/report and stdlib metadata are used. The disposable
virtualenv contains distribution metadata, not installed packages: pip check
can validate the proposed graph without changing the developer's environment.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import venv
from email.message import Message
from importlib import metadata
from pathlib import Path
from tempfile import TemporaryDirectory


def _name(value: str) -> str:
    return re.sub(r"[-_.]+", "-", value).lower()


def _package_metadata(name: str, version: str, requirements: list[str]) -> Message:
    message = Message()
    message["Metadata-Version"] = "2.1"
    message["Name"] = name
    message["Version"] = version
    for requirement in requirements:
        message["Requires-Dist"] = requirement
    return message


def _projected_metadata(report: dict) -> dict[str, Message]:
    if report["version"] != "1":
        raise ValueError(f"Unsupported pip installation report version: {report['version']}")
    proposed = {}
    for dist in metadata.distributions():
        # Match importlib's first-distribution precedence if paths overlap.
        name = dist.metadata["Name"]
        proposed.setdefault(_name(name), _package_metadata(name, dist.version, dist.requires or []))
    for item in report["install"]:
        fields = item["metadata"]
        proposed[_name(fields["name"])] = _package_metadata(
            fields["name"], fields["version"], fields.get("requires_dist", [])
        )
    return proposed


def _check_projection(proposed: dict[str, Message], directory: Path) -> None:
    environment = directory / "check-env"
    venv.EnvBuilder(with_pip=False, symlinks=sys.platform != "win32").create(environment)
    interpreter = environment / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    # A developer's PYTHONPATH must not shadow the proposed metadata with
    # distributions from the real environment and make an invalid plan pass.
    isolated = {
        key: value for key, value in os.environ.items() if key not in ("PYTHONPATH", "PYTHONHOME")
    }
    result = subprocess.run(
        [str(interpreter), "-c", "import sysconfig; print(sysconfig.get_path('purelib'))"],
        check=True,
        capture_output=True,
        text=True,
        env=isolated,
    )
    site_packages = Path(result.stdout.strip())
    for name, message in proposed.items():
        dist_info = site_packages / f"{name.replace('-', '_')}-{message['Version']}.dist-info"
        dist_info.mkdir()
        (dist_info / "METADATA").write_text(message.as_string(), encoding="utf-8")
    subprocess.run(
        [sys.executable, "-m", "pip", "--python", str(environment), "check"],
        check=True,
        env=isolated,
    )


def main(arguments: list[str]) -> int:
    try:
        pip_version = metadata.version("pip")
    except metadata.PackageNotFoundError:
        pip_version = "0"
    if int(pip_version.split(".", 1)[0]) < 23:
        print(
            "Development installs require pip >=23.0. "
            "Run `python -m ensurepip --upgrade` or `python -m pip install --upgrade pip`.",
            file=sys.stderr,
        )
        return 1

    install = [sys.executable, "-m", "pip", "install", *arguments]
    with TemporaryDirectory(prefix="tsarina-develop-") as temporary:
        directory = Path(temporary)
        report_path = directory / "report.json"
        print("Resolving all editable checkouts before installing ...", flush=True)
        try:
            subprocess.run([*install, "--dry-run", "--report", str(report_path)], check=True)
            proposed = _projected_metadata(json.loads(report_path.read_text(encoding="utf-8")))
            print("Checking proposed versions against installed consumers ...", flush=True)
            _check_projection(proposed, directory)
        except (subprocess.CalledProcessError, ValueError) as error:
            if isinstance(error, ValueError):
                print(str(error), file=sys.stderr)
            print(
                "Development install refused; installed packages were not changed. "
                "Update the incompatible sibling checkout or resolve the requirements shown above; "
                "use a separate virtualenv if consumers need different versions.",
                file=sys.stderr,
            )
            return error.returncode if isinstance(error, subprocess.CalledProcessError) else 1

        # Keep the second resolution at the versions checked above, including
        # dependencies already satisfied in the current environment.
        constraints = directory / "constraints.txt"
        constraints.write_text(
            "".join(
                f"{name}=={message['Version']}\n" for name, message in sorted(proposed.items())
            ),
            encoding="utf-8",
        )
        try:
            subprocess.run([*install, "--constraint", str(constraints)], check=True)
            print("Checking installed dependencies ...", flush=True)
            subprocess.run([sys.executable, "-m", "pip", "check"], check=True)
        except subprocess.CalledProcessError as error:
            print(
                "Development install failed; resolve the errors above before using it.",
                file=sys.stderr,
            )
            return error.returncode
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
