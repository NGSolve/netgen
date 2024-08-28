from subprocess import check_output
import os
import platform
import requests
import sys
import argparse


def is_package_available(package_name, version):
    architecture = platform.machine()
    py_version = "cp" + "".join(platform.python_version_tuple()[:2])

    url = f"https://pypi.org/pypi/{package_name}/{version}/json"

    try:
        response = requests.get(url)
        if response.status_code != 200:
            return False

        data = response.json()

        for file_info in data["urls"]:
            name = file_info.get("filename", "")
            if name.endswith(".whl") and py_version in name and architecture in name:
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


def get_version(cwd):
    git_version = (
        check_output(["git", "describe", "--tags"], cwd=cwd).decode("utf-8").strip()
    )

    version = git_version[1:].split("-")
    if len(version) > 2:
        version = version[:2]
    if len(version) > 1:
        version = ".post".join(version)
        if is_dev_build():
            version += ".dev0"
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

    version = get_version(args.dir)
    if args.get_version:
        print(version)
    elif args.check_pip:
        if is_package_available(args.package, version):
            print(f"{args.package}=={version} is already on pypi")
            sys.exit(1)
    else:
        print("no action")


if __name__ == "__main__":
    main()
