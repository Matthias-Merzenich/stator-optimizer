import os
import shutil
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version


def parse_version(version_string):
    return tuple(map(int, (version_string.split('.') + ['0', '0'])[:3]))


def validate_lifelib(need_latest_lifelib=False):
    min_version = "2.5.9"
    repo_url = "https://gitlab.com/apgoucher/lifelib.git"
    clone_required = False
    try:
        current_version = version("python-lifelib")
        if (
            parse_version(current_version) < parse_version(min_version)
            and not os.path.exists("lifelib")
            and need_latest_lifelib
        ):
            print("Incompatible lifelib version: "
                  f"{current_version} < {min_version}", file=sys.stderr
                 )
            clone_required = True
    except PackageNotFoundError:
        if not os.path.exists("lifelib"):
            print("lifelib not found.", file=sys.stderr)
            clone_required = True

    if clone_required:
        print("Fetching latest version...", file=sys.stderr)
        if shutil.which("git") is None:
                print("ERROR: 'git' is not installed or is not in PATH.",
                      file=sys.stderr)
                print("Please install git or download the "
                      "latest version of lifelib from\n"
                      "https://gitlab.com/apgoucher/lifelib", file=sys.stderr
                     )
                sys.exit(1)
        try:
            subprocess.run(["git", "clone", repo_url, "lifelib"], check=True)
            print("Clone successful.", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error cloning repository: {e}", file=sys.stderr)
            sys.exit(1)

    if os.path.exists("lifelib"):
        sys.path.insert(0, os.path.abspath("lifelib"))

    import lifelib
    return lifelib
