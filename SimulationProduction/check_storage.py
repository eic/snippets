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


def get_file_pfns(did, rse, max_files=None):
    """Get XRootD file paths for a DID at a specific RSE.

    Args:
        did: Rucio DID (can be dataset or file)
        rse: RSE name (BNL-XRD or EIC-XRD)
        max_files: Maximum number of files to query (None = all files)

    Returns:
        Tuple of (list of XRootD paths, total file count)
    """
    # XRootD prefix mapping for each RSE
    XROOTD_PREFIXES = {
        'BNL-XRD': 'root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC/',
        'EIC-XRD': 'root://dtn-eic.jlab.org//volatile/eic/EPIC'
    }

    try:
        # Get list of files in the dataset
        result = subprocess.run(
            ['/opt/local/bin/rucio', 'list-content', did],
            capture_output=True,
            text=True,
            timeout=60
        )

        # Extract file DIDs
        file_dids = []
        for line in result.stdout.strip().split('\n'):
            if 'FILE' in line:
                parts = line.split()
                if len(parts) >= 2:
                    file_dids.append(parts[1])

        total_files = len(file_dids)
        if not file_dids:
            return [], 0

        # Limit files if requested
        if max_files:
            file_dids = file_dids[:max_files]

        # Construct XRootD paths by replacing epic: prefix
        xrootd_prefix = XROOTD_PREFIXES.get(rse)
        if not xrootd_prefix:
            print(f"Error: Unknown RSE {rse}", file=sys.stderr)
            return [], 0

        pfns = []
        for file_did in file_dids:
            if file_did.startswith('epic:'):
                # Replace epic: with the XRootD prefix
                xrootd_path = xrootd_prefix + file_did[5:]  # Remove 'epic:'
                pfns.append(xrootd_path)

        return pfns, total_files

    except Exception as e:
        print(f"Error getting replicas for {did} at {rse}: {e}", file=sys.stderr)
        return [], 0


def main():
    parser = argparse.ArgumentParser(description='Check storage locations for rucio datasets')
    parser.add_argument('did', help='Rucio DID to check')
    parser.add_argument('--show-files', action='store_true', help='Show XRootD file paths')
    parser.add_argument('--rse', choices=['BNL-XRD', 'EIC-XRD'], help='Show files only for this RSE')
    parser.add_argument('--max-files', type=int, metavar='N', help='Show only first N files (default: 10)')
    parser.add_argument('--all-files', action='store_true', help='Show all files (can be slow for large datasets)')
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

        # Determine file limit
        if args.all_files:
            max_files = None
        elif args.max_files:
            max_files = args.max_files
        else:
            max_files = 10

        for rse in rses_to_check:
            print(f"\n{rse} files:")
            pfns, total = get_file_pfns(did, rse, max_files)
            if pfns:
                for pfn in pfns:
                    print(f"  {pfn}")
                if max_files and total > max_files:
                    print(f"  ... ({total - max_files} more files not shown, use --max-files or --all-files)")
            else:
                print(f"  (no files found)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
