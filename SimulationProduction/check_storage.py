#!/usr/bin/env python3
"""Check which storage locations (BNL-XRD or EIC-XRD) contain datasets."""

import subprocess
import sys
import argparse


def get_storage_locations(did):
    """Get storage locations for a rucio DID.

    Args:
        did: Rucio DID like epic:/RECO/26.04.1/epic_craterlake/DIS/...

    Returns:
        List of storage locations (BNL-XRD, EIC-XRD), excluding JLAB-TAPE-SE
    """
    try:
        result = subprocess.run(
            ['/opt/local/bin/rucio', 'list-rules', did],
            capture_output=True,
            text=True,
            timeout=30
        )

        locations = []
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            # Skip header and empty lines
            if not line or line.startswith('ID') or line.startswith('---'):
                continue

            # Parse the RSE_EXPRESSION column (4th column)
            parts = line.split()
            if len(parts) >= 5:
                rse = parts[4]
                if rse in ['BNL-XRD', 'EIC-XRD']:
                    locations.append(rse)

        return locations

    except Exception as e:
        print(f"Error checking rules for {did}: {e}", file=sys.stderr)
        return []


def get_file_pfns(did, rse):
    """Get XRootD file paths for a DID at a specific RSE.

    Args:
        did: Rucio DID
        rse: RSE name (BNL-XRD or EIC-XRD)

    Returns:
        List of XRootD paths
    """
    try:
        result = subprocess.run(
            ['/opt/local/bin/rucio', 'replica', 'list', 'file', '--pfns', '--rses', rse, did],
            capture_output=True,
            text=True,
            timeout=60
        )

        pfns = []
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if line.startswith('root://'):
                pfns.append(line)

        return pfns

    except Exception as e:
        print(f"Error getting replicas for {did} at {rse}: {e}", file=sys.stderr)
        return []


def main():
    parser = argparse.ArgumentParser(description='Check storage locations for rucio datasets')
    parser.add_argument('did', help='Rucio DID to check')
    parser.add_argument('--show-files', action='store_true', help='Show XRootD file paths')
    parser.add_argument('--rse', choices=['BNL-XRD', 'EIC-XRD'], help='Show files only for this RSE')
    args = parser.parse_args()

    did = args.did

    print(f"Checking storage for: {did}")
    locations = get_storage_locations(did)

    if not locations:
        print("  ✗ Not found in BNL-XRD or EIC-XRD")
        return 1

    print(f"  ✓ Found in: {', '.join(locations)}")

    if args.show_files:
        rses_to_check = [args.rse] if args.rse else locations

        for rse in rses_to_check:
            print(f"\n{rse} files:")
            pfns = get_file_pfns(did, rse)
            if pfns:
                for pfn in pfns:
                    print(f"  {pfn}")
            else:
                print(f"  (no files found)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
