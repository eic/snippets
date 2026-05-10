#!/usr/bin/env python3
"""Check if datasets are available in a specific campaign via rucio."""

import csv
import subprocess
import sys
import argparse
import urllib.request
import io

def check_rucio_did(path, campaign):
    """Check if DID exists in specified campaign for given path.

    Args:
        path: Path like /volatile/eic/EPIC/EVGEN/DDIS/rapgap3.310-1.0/noRad/ep/10x100
        campaign: Campaign version like "26.04" or "26.03"

    Returns:
        Tuple of (availability_marker, list_of_dids)
        availability_marker: 'X' if found, '' if not found
        list_of_dids: List of DID strings
    """
    # Extract portion after /volatile/eic/EPIC/EVGEN
    if not path.startswith('/volatile/eic/EPIC/EVGEN'):
        return '', []

    suffix = path.replace('/volatile/eic/EPIC/EVGEN', '').lstrip('/')
    if not suffix:
        return '', []

    # Special handling for UPSILON_ABCONV - check at file level
    if '/EXCLUSIVE/UPSILON_ABCONV/' in path:
        return check_upsilon_abconv(path, campaign)

    # Construct rucio DID pattern
    did_pattern = f"epic:/RECO*{campaign}*/{suffix}"

    try:
        result = subprocess.run(
            ['/opt/local/bin/rucio', 'did', 'list', '--short', did_pattern],
            capture_output=True,
            text=True,
            timeout=30
        )

        # Extract DIDs from output
        dids = []
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if line and line.startswith('epic:/RECO/'):
                dids.append(line)

        if dids:
            return 'X', dids
        return '', []
    except Exception as e:
        print(f"Error checking {did_pattern}: {e}", file=sys.stderr)
        return '', []


def check_upsilon_abconv(path, campaign):
    """Special handling for UPSILON_ABCONV datasets at file level.

    Args:
        path: Path like /volatile/eic/EPIC/EVGEN/EXCLUSIVE/UPSILON_ABCONV/upsilon1s_threshold_ab_hiAcc_10x100.hepmc3.tree.root
        campaign: Campaign version like "26.04"

    Returns:
        Tuple of (availability_marker, list_of_dids)
    """
    import os

    # Extract filename from path
    filename = os.path.basename(path)
    # Remove .hepmc3.tree.root extension to get base name
    if not filename.endswith('.hepmc3.tree.root'):
        return '', []

    base_name = filename.replace('.hepmc3.tree.root', '')

    # Find container DIDs for this campaign
    container_pattern = f"epic:*{campaign}*/EXCLUSIVE/UPSILON_ABCONV"

    try:
        result = subprocess.run(
            ['/opt/local/bin/rucio', 'did', 'list', '--short', container_pattern],
            capture_output=True,
            text=True,
            timeout=30
        )

        containers = []
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if line and line.startswith('epic:/RECO/'):
                containers.append(line)

        if not containers:
            return '', []

        # Check file contents of each container
        matching_files = []
        for container in containers:
            try:
                content_result = subprocess.run(
                    ['/opt/local/bin/rucio', 'did', 'content', 'list', '--short', container],
                    capture_output=True,
                    text=True,
                    timeout=30
                )

                for file_line in content_result.stdout.strip().split('\n'):
                    file_line = file_line.strip()
                    # Check if this file matches our base name
                    if base_name in file_line and file_line.endswith('.eicrecon.edm4eic.root'):
                        matching_files.append(file_line)

            except Exception as e:
                print(f"Error checking container {container}: {e}", file=sys.stderr)
                continue

        if matching_files:
            return 'X', matching_files
        return '', []

    except Exception as e:
        print(f"Error checking UPSILON_ABCONV pattern {container_pattern}: {e}", file=sys.stderr)
        return '', []

def main():
    parser = argparse.ArgumentParser(description='Check dataset availability in specific campaign')
    parser.add_argument('--campaign', default='26.04', help='Campaign version (e.g., 26.04, 26.03)')
    parser.add_argument('--output', default='datasets_with_campaign.csv', help='Output CSV file path')
    args = parser.parse_args()

    campaign = args.campaign
    column_name = f'In {campaign}?'

    output_file = args.output
    temp_file = output_file + '.tmp'

    # Fetch CSV from GitHub
    url = 'https://raw.githubusercontent.com/eic/epic-prod/main/docs/_data/datasets.csv'
    print(f"Fetching {url}")
    with urllib.request.urlopen(url) as response:
        csv_content = response.read().decode('utf-8')

    # Parse CSV from string
    csv_file = io.StringIO(csv_content)
    reader = csv.DictReader(csv_file)
    fieldnames = reader.fieldnames

    # Insert column after Priority column (before Dataset Path)
    priority_index = fieldnames.index('Priority')
    new_fieldnames = list(fieldnames[:priority_index+1]) + [column_name] + list(fieldnames[priority_index+1:])

    rows = list(reader)

    print(f"Checking {len(rows)} datasets for campaign {campaign}...")

    # Check each dataset
    for i, row in enumerate(rows, 1):
        path = row.get('Dataset Path', '')
        print(f"[{i}/{len(rows)}] Checking {path}...", file=sys.stderr)

        availability, dids = check_rucio_did(path, campaign)
        row[column_name] = availability

        if availability:
            print(f"  ✓ Found in {campaign}:", file=sys.stderr)
            for did in dids:
                print(f"    - {did}", file=sys.stderr)
        else:
            print(f"  ✗ Not in {campaign}", file=sys.stderr)

    # Write output
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=new_fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote results to {output_file}")

if __name__ == '__main__':
    main()
