#!/usr/bin/env python3
"""
Visualise the release timeline of a software component and its integration
into the ePIC stack container.

Reads the three output files produced by versions.py:
  - tags: all upstream release tags with dates
  - eic_container_tags: filtered eic/containers stable tags with dates
  - result: upstream version tag + date for each matched container tag

Usage:
  python plot.py [--component NAME] [--skip N] [--input-dir DIR] [--output FILE]
"""

import argparse
import os
from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


# Default number of early upstream tags to skip per component (to focus on
# the era relevant to ePIC/EIC containers).
DEFAULT_SKIP = {
    "Acts": 101,
    "DD4hep": 61,
}

# Annotation offset parameters per component.  Each callable receives the
# (integer) tag index and returns an (x, y) offset in points.
DEFAULT_ANNOTATION_PARAMS: dict[str, dict] = {
    "Acts": {
        "divisor": 60.0,
        "label_xy": lambda ix: ((-200 + ix * 1.5), (ix * 0.8 - 50)),
        "container_xy": lambda ix: ((-100 + ix), (-80 + ix * 0.5)),
    },
    "DD4hep": {
        "divisor": 8.0,
        "label_xy": lambda ix: ((-100 + 2 * ix), (ix * 2.4 - 10)),
        "container_xy": lambda ix: ((-10 + 2 * ix), (-80 + ix)),
    },
}


@np.vectorize
def to_datetime(date: str) -> datetime:
    return datetime.strptime(date, "%b %d %Y %H:%M:%S")


def annotate(ax, label, x, y, xytext):
    ax.annotate(
        label,
        xy=(x, y),
        xytext=xytext,
        textcoords="offset points",
        fontsize=5,
        arrowprops={"arrowstyle": "-|>", "color": "black", "linewidth": 1.0},
    )


def load_tag_file(path: str) -> np.ndarray:
    """Read a tag file (tag|date lines) and return a 2-column string array."""
    with open(path) as fp:
        rows = [l.strip().split("|") for l in fp if l.strip()]
    return np.array(rows)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--component", "-c",
        default="DD4hep",
        help="Component display name for axis labels (default: DD4hep). "
             f"Auto-skip configured for: {', '.join(sorted(DEFAULT_SKIP))}.",
    )
    parser.add_argument(
        "--skip", "-s",
        type=int,
        default=None,
        metavar="N",
        help="Number of early upstream tags to skip before plotting "
             "(default: auto per --component, 0 for unknown components).",
    )
    parser.add_argument(
        "--input-dir", "-i",
        default=".",
        metavar="DIR",
        help="Directory containing tags, eic_container_tags, result files "
             "(default: current directory).",
    )
    parser.add_argument(
        "--output", "-O",
        default="fig.png",
        metavar="FILE",
        help="Output figure path (default: fig.png).",
    )
    args = parser.parse_args()

    component = args.component
    skip = args.skip if args.skip is not None else DEFAULT_SKIP.get(component, 0)

    def inp(filename):
        return os.path.join(args.input_dir, filename)

    tags = load_tag_file(inp("tags"))
    tags_dates = to_datetime(tags[:, 1])
    order = np.argsort(tags_dates)
    tags = tags[order][skip:]
    tags_dates = tags_dates[order][skip:]

    eic_container_tags = load_tag_file(inp("eic_container_tags"))
    eic_container_tags_dates = to_datetime(eic_container_tags[:, 1])
    order = np.argsort(eic_container_tags_dates)
    eic_container_tags = eic_container_tags[order]
    eic_container_tags_dates = eic_container_tags_dates[order]

    result = load_tag_file(inp("result"))
    result_dates = to_datetime(result[:, 1])
    order = np.argsort(result_dates)
    result = result[order]
    result_dates = result_dates[order]

    if len(eic_container_tags) != len(result):
        raise ValueError(
            f"Length mismatch: {len(eic_container_tags)} entries in 'eic_container_tags' "
            f"but {len(result)} entries in 'result'. "
            "Re-run versions.py to regenerate consistent data files."
        )

    result_ixs, tags_ixs = np.nonzero(
        result[:, 0][:, np.newaxis] == tags[:, 0][np.newaxis, :]
    )
    print(result_ixs)
    print(tags_ixs)

    fig, ax = plt.subplots(figsize=(6.4, 3), layout="constrained")
    ax.plot(tags_dates, np.arange(len(tags)), label=f"{component} releases")
    ax.plot(
        eic_container_tags_dates[result_ixs],
        tags_ixs,
        label="ePIC stack releases",
    )

    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    ax.xaxis.set_major_locator(mpl.dates.YearLocator())
    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%Y"))

    container_bump_ix = np.nonzero(result[1:, 0] != result[:-1, 0])[0] + 1
    assert np.allclose(result_ixs, range(len(result_ixs)))

    ann_params = DEFAULT_ANNOTATION_PARAMS.get(component)
    for cix, ix in zip(container_bump_ix, tags_ixs[container_bump_ix]):
        if ann_params is None:
            continue
        f = np.tanh(ix / ann_params["divisor"])
        lxy = ann_params["label_xy"](ix)
        cxy = ann_params["container_xy"](ix)
        annotate(ax, tags[ix, 0], tags_dates[ix], ix,
                 (lxy[0] * f, lxy[1] * f))
        annotate(ax, eic_container_tags[cix, 0], eic_container_tags_dates[cix], ix,
                 (cxy[0] * f, cxy[1] * f))

    ax.set_ylabel(f"{component} version index", loc="top")
    ax.legend()
    ax.spines[["right", "top"]].set_visible(False)

    plt.savefig(args.output, dpi=300)
    print(f"Figure saved to {args.output}")
    plt.show()


if __name__ == "__main__":
    main()
