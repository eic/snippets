#!/usr/bin/env python3

"""
merge_stl.py

Perform a boolean union of two overlapping STL meshes and write
the result to a new STL file.

Usage:
    python merge_stl.py mesh1.stl mesh2.stl merged.stl

Optional:
    --engine manifold
    --engine blender
    --engine scad
"""

import argparse
import sys
from pathlib import Path

import trimesh


def load_mesh(path):
    """Load an STL file as a Trimesh."""
    mesh = trimesh.load_mesh(
        path,
        force="mesh",
        process=True
    )

    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError(
            "{} did not load as a single mesh.".format(path)
        )

    return mesh


def combine_result(result):
    """
    Convert a boolean result into a single Trimesh.
    """
    if result is None:
        return None

    if isinstance(result, trimesh.Trimesh):
        return result

    if isinstance(result, trimesh.Scene):
        meshes = [
            geom
            for geom in result.geometry.values()
            if isinstance(geom, trimesh.Trimesh)
        ]

        if not meshes:
            return None

        return trimesh.util.concatenate(meshes)

    if isinstance(result, (list, tuple)):
        meshes = [
            geom
            for geom in result
            if isinstance(geom, trimesh.Trimesh)
        ]

        if not meshes:
            return None

        return trimesh.util.concatenate(meshes)

    raise TypeError(
        "Unsupported boolean result type: {}".format(type(result))
    )


def merge_meshes(mesh_a, mesh_b, engine=None):
    """Perform boolean union of two meshes."""
    result = trimesh.boolean.union(
        [mesh_a, mesh_b],
        engine=engine
    )

    return combine_result(result)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge two STL meshes using a boolean union."
        )
    )

    parser.add_argument(
        "input_a",
        help="First STL mesh"
    )

    parser.add_argument(
        "input_b",
        help="Second STL mesh"
    )

    parser.add_argument(
        "output",
        help="Output STL mesh"
    )

    parser.add_argument(
        "--engine",
        default=None,
        help=(
            "Boolean engine: manifold, blender, or scad. "
            "Default: auto"
        )
    )

    args = parser.parse_args()

    # Load meshes
    try:
        mesh_a = load_mesh(args.input_a)
        mesh_b = load_mesh(args.input_b)
    except Exception as exc:
        print(
            "Error loading mesh: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(1)

    print("Mesh A:")
    print("  vertices:   {}".format(len(mesh_a.vertices)))
    print("  faces:      {}".format(len(mesh_a.faces)))
    print("  watertight: {}".format(mesh_a.is_watertight))
    print("  volume:     {}".format(abs(mesh_a.volume)))

    print()
    print("Mesh B:")
    print("  vertices:   {}".format(len(mesh_b.vertices)))
    print("  faces:      {}".format(len(mesh_b.faces)))
    print("  watertight: {}".format(mesh_b.is_watertight))
    print("  volume:     {}".format(abs(mesh_b.volume)))

    if not mesh_a.is_watertight or not mesh_b.is_watertight:
        print()
        print(
            "WARNING: one or both meshes are not watertight. "
            "Boolean union may fail or produce invalid geometry.",
            file=sys.stderr
        )

    # Boolean union
    print()
    print("Calculating boolean union...")

    try:
        result = merge_meshes(
            mesh_a,
            mesh_b,
            engine=args.engine
        )
    except Exception as exc:
        print(
            "Boolean union failed: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(2)

    if result is None or len(result.faces) == 0:
        print(
            "Boolean union produced an empty mesh.",
            file=sys.stderr
        )
        sys.exit(2)

    # Basic cleanup
#    result.remove_duplicate_faces()
#    result.remove_degenerate_faces()
    result.update_faces(result.unique_faces())
    result.remove_unreferenced_vertices()

    # Write output
    output_path = Path(args.output)

    try:
        result.export(str(output_path))
    except Exception as exc:
        print(
            "Error writing output: {}".format(exc),
            file=sys.stderr
        )
        sys.exit(1)

    print()
    print("Merged mesh:")
    print("  vertices:   {}".format(len(result.vertices)))
    print("  faces:      {}".format(len(result.faces)))
    print("  watertight: {}".format(result.is_watertight))
    print("  volume:     {}".format(abs(result.volume)))
    print()
    print("Wrote:")
    print("  {}".format(output_path))


if __name__ == "__main__":
    main()
