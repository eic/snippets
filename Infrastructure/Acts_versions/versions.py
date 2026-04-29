#!/usr/bin/env python3
"""
Extract package versions used in eic/containers repository tags.

Instead of running singularity/spack inside containers, this script checks out
the eic/containers repo at each tag and reads the package version directly from
either spack-environment/packages.yaml (newer format) or
spack-environment/dev/spack.yaml (older format) or spack.yaml (earliest tags).

Produces the same output files as versions.sh / versions_dd4hep.sh:
  - tags: all upstream release tags with dates
  - eic_container_tags: filtered eic/containers tags with dates
  - result: upstream version tag + date for each container tag

Usage:
  python versions.py [package]

  package: acts (default), dd4hep, or any spack package name
"""

import argparse
import os
import re
import subprocess
import sys

import yaml


# ---------------------------------------------------------------------------
# Package configurations
#
# Each entry maps a spack package name to:
#   repo_url:           upstream git repository
#   repo_dir:           local clone directory name
#   version_to_tag:     convert spack version string to upstream git tag
# ---------------------------------------------------------------------------

def _identity_version_to_tag(version):
    """v{version} — used by acts and most projects."""
    return f"v{version}"


def _dd4hep_version_to_tag(version):
    """
    DD4hep spack version '1.25.1' -> git tag 'v01-25-01'.
    Each dot-separated component is zero-padded to 2 digits, joined by '-'.
    """
    parts = version.split(".")
    return "v" + "-".join(p.zfill(2) for p in parts)


PACKAGES = {
    "acts": {
        "repo_url": "https://github.com/acts-project/acts",
        "repo_dir": "acts",
        "version_to_tag": _identity_version_to_tag,
    },
    "dd4hep": {
        "repo_url": "https://github.com/AIDASoft/DD4hep.git",
        "repo_dir": "DD4hep",
        "version_to_tag": _dd4hep_version_to_tag,
    },
}


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------

def run(cmd, **kwargs):
    """Run a shell command and return stdout."""
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        print(f"Command failed: {' '.join(cmd)}", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
    return result.stdout


def clone_if_missing(url, path):
    """Clone a git repository if the directory doesn't already exist."""
    if not os.path.isdir(path):
        print(f"Cloning {url}...")
        run(["git", "clone", url, path])
    else:
        print(f"Using existing {path}")


def get_tags_with_dates(repo_path):
    """
    Get all tags with their creator dates from a git repo.
    Returns list of "tag|date" strings.
    """
    output = run(
        ["git", "-C", repo_path, "for-each-ref",
         "--format=%(refname:short)|%(creatordate:format:%b %d %Y %H:%M:%S)",
         "refs/tags/*"]
    )
    return [line for line in output.strip().splitlines() if line]


def git_show(repo_path, ref, filepath):
    """Show file content at a given git ref. Returns None if file doesn't exist."""
    result = subprocess.run(
        ["git", "-C", repo_path, "show", f"{ref}:{filepath}"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return None
    return result.stdout


# ---------------------------------------------------------------------------
# Version extraction from container files
# ---------------------------------------------------------------------------

def extract_version_from_spack_yaml(content, spack_name):
    """
    Extract version from spack.yaml / spack-environment/dev/spack.yaml format:
      spack:
        specs:
          - <spack_name>@<version> ...
    """
    data = yaml.safe_load(content)
    if not data:
        return None
    specs = data.get("spack", {}).get("specs", [])
    # Match exact package name followed by @version (avoid e.g. acts-dd4hep matching dd4hep)
    pattern = re.compile(r"^" + re.escape(spack_name) + r"@([^\s]+)")
    for spec in specs:
        m = pattern.match(spec)
        if m:
            return m.group(1)
    return None


def extract_version_from_packages_yaml(content, spack_name):
    """
    Extract version from packages.yaml format:
      packages:
        <spack_name>:
          require:
          - '@<version>'
          - ...
    """
    data = yaml.safe_load(content)
    if not data:
        return None
    pkg = data.get("packages", {}).get(spack_name, {})
    require = pkg.get("require", [])
    for entry in require:
        if isinstance(entry, str):
            m = re.match(r"'?@([0-9][^']*)'?", entry)
            if m:
                return m.group(1)
        elif isinstance(entry, dict):
            spec = entry.get("spec", "")
            m = re.match(r"'?@([0-9][^']*)'?", spec)
            if m:
                return m.group(1)
    return None


def get_package_version(containers_path, tag, spack_name):
    """
    Try to extract a spack package version from a container tag.
    Searches multiple file locations in order.
    """
    # Try packages.yaml first (newer format, ~v24.03+)
    content = git_show(containers_path, tag, "spack-environment/packages.yaml")
    if content:
        version = extract_version_from_packages_yaml(content, spack_name)
        if version:
            return version

    # Try dev/spack.yaml (middle era, ~v23.05–v24.02)
    content = git_show(containers_path, tag, "spack-environment/dev/spack.yaml")
    if content:
        version = extract_version_from_spack_yaml(content, spack_name)
        if version:
            return version

    # Try spack.yaml or spack.yml at root (earliest tags, ~v23.03)
    for filename in ("spack.yaml", "spack.yml"):
        content = git_show(containers_path, tag, filename)
        if content:
            version = extract_version_from_spack_yaml(content, spack_name)
            if version:
                return version

    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract package versions from eic/containers tags."
    )
    parser.add_argument(
        "package",
        nargs="?",
        default="acts",
        help="Spack package name (default: acts). "
             f"Preconfigured: {', '.join(sorted(PACKAGES))}. "
             "Other names will be looked up with v{{version}} tag convention.",
    )
    args = parser.parse_args()

    spack_name = args.package
    pkg_config = PACKAGES.get(spack_name, {
        "repo_url": None,
        "repo_dir": None,
        "version_to_tag": _identity_version_to_tag,
    })

    repo_url = pkg_config["repo_url"]
    repo_dir = pkg_config["repo_dir"]
    version_to_tag = pkg_config["version_to_tag"]

    if repo_url is None:
        print(f"No upstream repository configured for '{spack_name}'.", file=sys.stderr)
        print(f"Known packages: {', '.join(sorted(PACKAGES))}", file=sys.stderr)
        print("Add an entry to PACKAGES in this script, or specify a known package.", file=sys.stderr)
        sys.exit(1)

    # Clone repos
    clone_if_missing(repo_url, repo_dir)
    clone_if_missing("https://github.com/eic/containers.git", "containers")

    # Get all upstream tags with dates
    print(f"Fetching {spack_name} tags from {repo_dir}...")
    upstream_tags_lines = get_tags_with_dates(repo_dir)
    with open("tags", "w") as f:
        for line in upstream_tags_lines:
            print(line)
            f.write(line + "\n")

    # Build lookup: tag_name -> "tag|date" line
    upstream_tag_lookup = {}
    for line in upstream_tags_lines:
        tag, date = line.split("|", 1)
        upstream_tag_lookup[tag] = line

    # Get container tags (all stable tags)
    print("Fetching container tags...")
    all_container_lines = get_tags_with_dates("containers")
    eic_container_tags = []
    for line in all_container_lines:
        tag = line.split("|")[0]
        if "stable" not in tag:
            continue
        eic_container_tags.append(line)

    with open("eic_container_tags", "w") as f:
        for line in eic_container_tags:
            print(line)
            f.write(line + "\n")

    # For each container tag, extract the package version and look up its upstream tag date
    print(f"Extracting {spack_name} versions from container tags...")
    results = []
    for line in eic_container_tags:
        container_tag = line.split("|")[0]
        version = get_package_version("containers", container_tag, spack_name)
        if version is None:
            print(f"  {container_tag}: {spack_name} version not found", file=sys.stderr)
            continue

        upstream_tag = version_to_tag(version)
        if upstream_tag in upstream_tag_lookup:
            result_line = upstream_tag_lookup[upstream_tag]
            results.append(result_line)
            print(f"  {container_tag} -> {spack_name} {version}: {result_line}")
        else:
            print(f"  {container_tag} -> {spack_name} {version}: "
                  f"tag {upstream_tag} not found in {repo_dir} repo",
                  file=sys.stderr)

    with open("result", "w") as f:
        for line in results:
            print(line)
            f.write(line + "\n")

    print(f"\nDone. {len(results)} results written.")


if __name__ == "__main__":
    main()
