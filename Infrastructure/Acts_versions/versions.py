#!/usr/bin/env python3
"""
Extract package versions used in eic/containers repository tags.

Instead of running singularity/spack inside containers, this script checks out
the eic/containers repo at each tag and reads the package version directly from
either spack-environment/packages.yaml (newer format) or
spack-environment/dev/spack.yaml (older format) or spack.yaml (earliest tags).

Unlike versions.sh / versions_dd4hep.sh (which require singularity and skip
early container tags), this script works purely from git history and covers
all stable tags including v23.* and v24.02.*.

Output files:
  - tags: all upstream release tags with dates
  - eic_container_tags: filtered eic/containers stable tags with dates
  - result: upstream version tag + date for each matched container tag

Usage:
  python versions.py [package] [--output-dir DIR] [--no-fetch]
  python versions.py --list-packages

  package: acts (default), dd4hep, eicrecon, podio, edm4hep, geant4,
           or any other preconfigured package name.
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
    """v{version} — used by acts, eicrecon, and most semver projects."""
    return f"v{version}"


def _hep_version_to_tag(version):
    """
    HEP-style spack version '1.25.1' -> git tag 'v01-25-01'.
    Each dot-separated component is zero-padded to 2 digits, joined by '-'.
    Used by DD4hep, podio, and EDM4hep.
    """
    parts = version.split(".")
    return "v" + "-".join(p.zfill(2) for p in parts)


def _geant4_version_to_tag(version):
    """
    Geant4 spack version '11.3.2' or '11.3.2.east' -> git tag 'v11.3.2'.
    Strips EIC-specific patch suffixes (e.g., '.east') beyond the 3rd component.
    Note: versions ending in '.east' are EIC-patched builds with no exact upstream tag;
    the tag for the unpatched base release is returned instead.
    """
    parts = version.split(".")
    return "v" + ".".join(parts[:3])


PACKAGES = {
    "acts": {
        "repo_url": "https://github.com/acts-project/acts",
        "repo_dir": "acts",
        "version_to_tag": _identity_version_to_tag,
    },
    "dd4hep": {
        "repo_url": "https://github.com/AIDASoft/DD4hep",
        "repo_dir": "DD4hep",
        "version_to_tag": _hep_version_to_tag,
    },
    "eicrecon": {
        "repo_url": "https://github.com/eic/EICrecon",
        "repo_dir": "EICrecon",
        "version_to_tag": _identity_version_to_tag,
    },
    "podio": {
        "repo_url": "https://github.com/AIDASoft/podio",
        "repo_dir": "podio",
        "version_to_tag": _hep_version_to_tag,
    },
    "edm4hep": {
        "repo_url": "https://github.com/key4hep/EDM4hep",
        "repo_dir": "EDM4hep",
        "version_to_tag": _hep_version_to_tag,
    },
    "geant4": {
        "repo_url": "https://github.com/Geant4/geant4",
        "repo_dir": "geant4",
        "version_to_tag": _geant4_version_to_tag,
    },
}


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------

def run(cmd, **kwargs):
    """Run a shell command and return the CompletedProcess result."""
    return subprocess.run(cmd, capture_output=True, text=True, **kwargs)


def clone_or_fetch(url, path, fetch=True):
    """
    Clone a git repository if the directory doesn't exist; otherwise fetch new tags.
    Exits with an error if the initial clone fails.
    """
    if not os.path.isdir(path):
        print(f"Cloning {url} -> {path} ...")
        result = run(["git", "clone", url, path])
        if result.returncode != 0:
            print(f"ERROR: git clone failed for {url}:", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            sys.exit(1)
    elif fetch:
        print(f"Fetching new tags in {path} ...")
        result = run(["git", "-C", path, "fetch", "--tags", "--quiet"])
        if result.returncode != 0:
            print(
                f"WARNING: git fetch failed for {path}: {result.stderr.strip()}",
                file=sys.stderr,
            )


def get_tags_with_dates(repo_path):
    """
    Get all tags with their creator dates from a git repo.
    Returns list of "tag|date" strings.
    Exits with an error if the git command fails.
    """
    result = run(
        ["git", "-C", repo_path, "for-each-ref",
         "--format=%(refname:short)|%(creatordate:format:%b %d %Y %H:%M:%S)",
         "refs/tags/*"],
    )
    if result.returncode != 0:
        print(f"ERROR: failed to list tags in {repo_path}:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(1)
    return [line for line in result.stdout.strip().splitlines() if line]


def git_show(repo_path, ref, filepath):
    """Show file content at a given git ref. Returns None if file doesn't exist at that ref."""
    result = run(["git", "-C", repo_path, "show", f"{ref}:{filepath}"])
    if result.returncode != 0:
        return None
    return result.stdout


# ---------------------------------------------------------------------------
# Version extraction from container files
# ---------------------------------------------------------------------------

def extract_version_from_spack_yaml(content, spack_name):
    """
    Extract version from spack.yaml / spack-environment/dev/spack.yaml format::

        spack:
          specs:
            - <spack_name>@<version> ...

    Uses an exact package-name match to avoid e.g. ``acts-dd4hep`` matching ``dd4hep``.
    """
    try:
        data = yaml.safe_load(content)
    except yaml.YAMLError:
        return None
    if not data:
        return None
    specs = data.get("spack", {}).get("specs", [])
    pattern = re.compile(r"^" + re.escape(spack_name) + r"@([^\s]+)")
    for spec in specs:
        if not isinstance(spec, str):
            continue
        m = pattern.match(spec)
        if m:
            return m.group(1)
    return None


def extract_version_from_packages_yaml(content, spack_name):
    """
    Extract version from packages.yaml format::

        packages:
          <spack_name>:
            require:
            - '@<version>'

    Handles both plain string entries and ``{spec: ..., when: ...}`` dict entries.
    After yaml.safe_load(), YAML quoting is already stripped, so the entry is
    a plain Python string such as ``@44.4.0``.
    """
    try:
        data = yaml.safe_load(content)
    except yaml.YAMLError:
        return None
    if not data:
        return None
    pkg = data.get("packages", {}).get(spack_name, {})
    if not pkg:
        return None
    require = pkg.get("require", [])
    version_re = re.compile(r"@([0-9][^\s]*)")
    for entry in require:
        if isinstance(entry, str):
            m = version_re.match(entry)
            if m:
                return m.group(1)
        elif isinstance(entry, dict):
            spec = entry.get("spec", "")
            m = version_re.match(spec)
            if m:
                return m.group(1)
    return None


def get_package_version(containers_path, tag, spack_name):
    """
    Try to extract a spack package version from a container tag.

    Searches the following paths in order (newest to oldest format):

    1. ``spack-environment/packages.yaml``  — v24.03+
    2. ``spack-environment/dev/spack.yaml`` — v23.05–v24.02
    3. ``spack.yaml``                       — earliest tags (~v23.03)
    4. ``spack.yml``                        — fallback spelling
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

    # Try root-level spack.yaml / spack.yml (earliest tags, ~v23.03)
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
        description="Extract package versions from eic/containers tags.",
    )
    parser.add_argument(
        "package",
        nargs="?",
        default="acts",
        help=(
            "Spack package name to track (default: acts). "
            f"Preconfigured: {', '.join(sorted(PACKAGES))}. "
            "Use --list-packages for details."
        ),
    )
    parser.add_argument(
        "--output-dir", "-o",
        default=".",
        metavar="DIR",
        help="Directory for output files (tags, eic_container_tags, result). Default: '.'",
    )
    parser.add_argument(
        "--no-fetch",
        action="store_true",
        help="Skip 'git fetch' when repos already exist (use cached data).",
    )
    parser.add_argument(
        "--list-packages",
        action="store_true",
        help="List all preconfigured packages and exit.",
    )
    args = parser.parse_args()

    if args.list_packages:
        print("Preconfigured packages:")
        for name, cfg in sorted(PACKAGES.items()):
            print(f"  {name:12s}  {cfg['repo_url']}")
        return

    spack_name = args.package
    if spack_name not in PACKAGES:
        print(f"ERROR: Unknown package '{spack_name}'.", file=sys.stderr)
        print(f"  Known packages: {', '.join(sorted(PACKAGES))}", file=sys.stderr)
        print("  Add an entry to PACKAGES in this script to support additional packages.", file=sys.stderr)
        sys.exit(1)

    pkg_config = PACKAGES[spack_name]
    repo_url = pkg_config["repo_url"]
    repo_dir = pkg_config["repo_dir"]
    version_to_tag = pkg_config["version_to_tag"]
    fetch = not args.no_fetch

    os.makedirs(args.output_dir, exist_ok=True)

    def out(filename):
        return os.path.join(args.output_dir, filename)

    # Clone / update repos
    clone_or_fetch(repo_url, repo_dir, fetch=fetch)
    clone_or_fetch("https://github.com/eic/containers", "containers", fetch=fetch)

    # Get all upstream tags with dates
    print(f"\nFetching {spack_name} tags from '{repo_dir}' ...")
    upstream_tags_lines = get_tags_with_dates(repo_dir)
    with open(out("tags"), "w") as f:
        for line in upstream_tags_lines:
            f.write(line + "\n")
    print(f"  {len(upstream_tags_lines)} upstream tags found.")

    # Build lookup: tag_name -> "tag|date" line
    upstream_tag_lookup: dict[str, str] = {}
    for line in upstream_tags_lines:
        tag, _ = line.split("|", 1)
        upstream_tag_lookup[tag] = line

    # Get container stable tags
    print("\nFetching eic/containers tags ...")
    all_container_lines = get_tags_with_dates("containers")
    eic_container_tags = [
        line for line in all_container_lines
        if "stable" in line.split("|")[0]
    ]
    with open(out("eic_container_tags"), "w") as f:
        for line in eic_container_tags:
            f.write(line + "\n")
    print(f"  {len(eic_container_tags)} stable container tags found.")

    # For each container tag, extract the package version and look up its upstream tag date
    print(f"\nExtracting {spack_name} versions from container tags ...")
    results = []
    missing = 0
    for line in eic_container_tags:
        container_tag = line.split("|")[0]
        version = get_package_version("containers", container_tag, spack_name)
        if version is None:
            print(
                f"  WARNING: {container_tag}: {spack_name} not found in spack environment",
                file=sys.stderr,
            )
            missing += 1
            continue

        upstream_tag = version_to_tag(version)
        if upstream_tag in upstream_tag_lookup:
            result_line = upstream_tag_lookup[upstream_tag]
            results.append(result_line)
            print(f"  {container_tag} -> {spack_name} {version} ({upstream_tag})")
        else:
            print(
                f"  WARNING: {container_tag} -> {spack_name} {version}: "
                f"tag '{upstream_tag}' not found in '{repo_dir}' repo",
                file=sys.stderr,
            )
            missing += 1

    with open(out("result"), "w") as f:
        for line in results:
            f.write(line + "\n")

    print(f"\nDone: {len(results)} results written to '{args.output_dir}/result'.")
    if missing:
        print(f"  {missing} container tag(s) had no matching upstream tag (see warnings above).")


if __name__ == "__main__":
    main()
