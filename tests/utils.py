import argparse
import os
import requests
import sys
import time
from subprocess import check_output
from packaging import tags
from packaging.utils import parse_wheel_filename


_sys_tags = None


def _is_wheel_compatible(wheel_filename: str):
    global _sys_tags
    try:
        if _sys_tags is None:
            _sys_tags = set(tags.sys_tags())

        for tag in parse_wheel_filename(wheel_filename)[-1]:
            if tag in _sys_tags:
                return True

        return False
    except Exception as e:
        print(f"Error parsing wheel file: {e}")
        return False


def is_package_available(package_name, version):
    url = f"https://pypi.org/pypi/{package_name}/{version}/json"

    try:
        response = requests.get(url)
        if response.status_code != 200:
            return False

        data = response.json()

        for file_info in data["urls"]:
            name = file_info.get("filename", "")
            if _is_wheel_compatible(name):
                return True

        return False

    except requests.RequestException as e:
        print(f"Error checking package: {e}")
        return False


def is_dev_build():
    if "NG_NO_DEV_PIP_VERSION" in os.environ:
        return False
    if (
        "CI_COMMIT_REF_NAME" in os.environ
        and os.environ["CI_COMMIT_REF_NAME"] == "release"
    ):
        return False
    return True


def get_git_version(cwd):
    return check_output(["git", "describe", "--tags"], cwd=cwd).decode("utf-8").strip()


def get_dev_extension(cwd):
    if not is_dev_build():
        return ""

    # if the current commit does not belong to master, build a .dev1 package to avoid name conflicts for subsequent nightly builds from master
    try:
        check_output(["git", "merge-base", "--is-ancestor", "HEAD", "master"], cwd=cwd)
        return ".dev0"
    except:
        return ".dev1"


def get_version(cwd):
    git_version = get_git_version(cwd)

    version = git_version[1:].split("-")
    if len(version) > 2:
        version = version[:2]
    if len(version) > 1:
        version = ".post".join(version)
        version += get_dev_extension(cwd)
    else:
        version = version[0]

    return version


def main():
    parser = argparse.ArgumentParser(description="Netgen pip building utilities")
    parser.add_argument(
        "--check-pip",
        action="store_true",
        help="Check if package is on pypi already, fails with exit code 1 if available",
    )
    parser.add_argument(
        "--wait-pip",
        action="store_true",
        help="Wait until package is on pypi, fails with exit code 1 if still not available after 300s",
    )
    parser.add_argument(
        "--get-git-version",
        action="store_true",
        help="Generate the current package git version string",
    )
    parser.add_argument(
        "--get-version",
        action="store_true",
        help="Generate the current package version using git",
    )
    parser.add_argument("--dir", type=str, default=".", help="CWD to run git commands")
    parser.add_argument(
        "--package",
        type=str,
        default="netgen-mesher",
        help="Package name to check on pypi",
    )

    args = parser.parse_args()

    if args.get_git_version:
        print(get_git_version(args.dir))
    elif args.get_version:
        print(get_version(args.dir))
    elif args.check_pip:
        version = get_version(args.dir)
        if is_package_available(args.package, version):
            print(f"{args.package}=={version} is already on pypi")
            sys.exit(1)
    elif args.wait_pip:
        version = get_version(args.dir)
        t0 = time.time()
        while time.time() - t0 < 300 and not is_package_available(
            args.package, version
        ):
            time.sleep(20)

        if not is_package_available(args.package, version):
            print(
                f"Timeout waiting for package {args.package}=={version} to be available on pypi"
            )
            sys.exit(1)
    else:
        print("no action")


if __name__ == "__main__":
    main()
