#!/usr/bin/env python3
"""Recover generator particles from an EDM4hep npsim output file.
Shujie Li with Codex, 2026.

The input file contains both the original generator record and particles
created later by Geant4 interactions in detector material.  This script keeps
only MCParticles for which BITCreatedInSimulation (simulatorStatus bit 30) is
not set, rebuilds their parent/daughter graph, and writes a HepMC3 ASCII file
that can be used as a new npsim input.

Dependencies: numpy, awkward, uproot, pyhepmc. No local modules or detector
geometry are needed. Input positions/time are mm/ns; momentum/mass are GeV.

Example:
    python extract_generator_hepmc.py \
        hepmc/edm4hep/sample.edm4hep.root \
        hepmc/sample.hepmc3
"""

import argparse
import math
from pathlib import Path
import sys

import awkward as ak
import numpy as np
import pyhepmc
import uproot


CREATED_IN_SIMULATION_BIT = np.uint32(1 << 30)
SPEED_OF_LIGHT_MM_PER_NS = 299.792458

BRANCHES = {
    "event_number": "EventHeader/EventHeader.eventNumber",
    "pdg": "MCParticles/MCParticles.PDG",
    "generator_status": "MCParticles/MCParticles.generatorStatus",
    "simulator_status": "MCParticles/MCParticles.simulatorStatus",
    "mass": "MCParticles/MCParticles.mass",
    "time": "MCParticles/MCParticles.time",
    "vertex_x": "MCParticles/MCParticles.vertex.x",
    "vertex_y": "MCParticles/MCParticles.vertex.y",
    "vertex_z": "MCParticles/MCParticles.vertex.z",
    "endpoint_x": "MCParticles/MCParticles.endpoint.x",
    "endpoint_y": "MCParticles/MCParticles.endpoint.y",
    "endpoint_z": "MCParticles/MCParticles.endpoint.z",
    "momentum_x": "MCParticles/MCParticles.momentum.x",
    "momentum_y": "MCParticles/MCParticles.momentum.y",
    "momentum_z": "MCParticles/MCParticles.momentum.z",
    "parents_begin": "MCParticles/MCParticles.parents_begin",
    "parents_end": "MCParticles/MCParticles.parents_end",
    "daughters_begin": "MCParticles/MCParticles.daughters_begin",
    "daughters_end": "MCParticles/MCParticles.daughters_end",
    "parent_index": "_MCParticles_parents/_MCParticles_parents.index",
    "daughter_index": "_MCParticles_daughters/_MCParticles_daughters.index",
}


def event_array(arrays, branch_name):
    """Return one event's jagged branch as a NumPy array."""
    return np.asarray(ak.to_numpy(arrays[branch_name][0]))


def source_relations(begin, end, relation_indices, source_indices, source_set):
    """Return retained relation indices for each retained particle."""
    relations = {}
    number_of_particles = len(begin)

    for particle_index in source_indices:
        relation_begin = int(begin[particle_index])
        relation_end = int(end[particle_index])
        if not 0 <= relation_begin <= relation_end <= len(relation_indices):
            raise ValueError(
                f"particle {particle_index} has an invalid relation range "
                f"[{relation_begin}, {relation_end})"
            )

        retained_relations = []
        for related_index in relation_indices[relation_begin:relation_end]:
            related_index = int(related_index)
            if not 0 <= related_index < number_of_particles:
                raise ValueError(
                    f"particle {particle_index} refers to invalid particle "
                    f"index {related_index}"
                )
            if related_index in source_set:
                if related_index == particle_index:
                    raise ValueError(f"particle {particle_index} refers to itself")
                retained_relations.append(related_index)

        if len(set(retained_relations)) != len(retained_relations):
            raise ValueError(f"particle {particle_index} has duplicate relations")
        # Relation order has no physical meaning. Canonical ordering ensures
        # children with the same parents are grouped into the same vertex.
        relations[int(particle_index)] = tuple(sorted(retained_relations))

    return relations


def build_hepmc_event(arrays, entry_number):
    """Convert one uproot event into a source-only HepMC3 GenEvent."""
    pdg = event_array(arrays, BRANCHES["pdg"])
    generator_status = event_array(arrays, BRANCHES["generator_status"])
    simulator_status = event_array(arrays, BRANCHES["simulator_status"])
    mass = event_array(arrays, BRANCHES["mass"])
    time = event_array(arrays, BRANCHES["time"])
    vertex_x = event_array(arrays, BRANCHES["vertex_x"])
    vertex_y = event_array(arrays, BRANCHES["vertex_y"])
    vertex_z = event_array(arrays, BRANCHES["vertex_z"])
    endpoint_x = event_array(arrays, BRANCHES["endpoint_x"])
    endpoint_y = event_array(arrays, BRANCHES["endpoint_y"])
    endpoint_z = event_array(arrays, BRANCHES["endpoint_z"])
    momentum_x = event_array(arrays, BRANCHES["momentum_x"])
    momentum_y = event_array(arrays, BRANCHES["momentum_y"])
    momentum_z = event_array(arrays, BRANCHES["momentum_z"])

    # Bit 30 marks particles made by detector simulation.  uint32 avoids
    # signed-integer surprises when testing the high status bit.
    source_mask = (
        simulator_status.astype(np.uint32) & CREATED_IN_SIMULATION_BIT
    ) == 0
    source_indices = np.flatnonzero(source_mask)
    source_set = set(source_indices.tolist())

    parents = source_relations(
        event_array(arrays, BRANCHES["parents_begin"]),
        event_array(arrays, BRANCHES["parents_end"]),
        event_array(arrays, BRANCHES["parent_index"]),
        source_indices,
        source_set,
    )
    daughters = source_relations(
        event_array(arrays, BRANCHES["daughters_begin"]),
        event_array(arrays, BRANCHES["daughters_end"]),
        event_array(arrays, BRANCHES["daughter_index"]),
        source_indices,
        source_set,
    )

    # A valid HepMC graph records each relation in both directions.  Stop on
    # inconsistent input rather than silently producing a different event.
    for child_index in source_indices:
        for parent_index in parents[int(child_index)]:
            if int(child_index) not in daughters[parent_index]:
                raise ValueError(
                    f"parent/daughter mismatch between particles "
                    f"{parent_index} and {child_index}"
                )
        for daughter_index in daughters[int(child_index)]:
            if int(child_index) not in parents[daughter_index]:
                raise ValueError(
                    f"parent/daughter mismatch between particles "
                    f"{child_index} and {daughter_index}"
                )

    event_number_array = event_array(arrays, BRANCHES["event_number"])
    event_number = (
        int(event_number_array[0]) if len(event_number_array) else entry_number
    )
    event = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
    event.event_number = event_number

    # Keep the original MCParticle order.  This makes direct comparisons to
    # the EDM4hep source subset straightforward.
    particles = {}
    for particle_index in source_indices:
        particle_index = int(particle_index)
        px = float(momentum_x[particle_index])
        py = float(momentum_y[particle_index])
        pz = float(momentum_z[particle_index])
        particle_mass = float(mass[particle_index])
        energy = math.sqrt(
            max(0.0, px * px + py * py + pz * pz + particle_mass * particle_mass)
        )

        particle = pyhepmc.GenParticle(
            pyhepmc.FourVector(px, py, pz, energy),
            int(pdg[particle_index]),
            int(generator_status[particle_index]),
        )
        particle.generated_mass = particle_mass
        event.add_particle(particle)
        particles[particle_index] = particle

    # Children with the same parent set share one HepMC interaction vertex.
    # Use the parent endpoint as the canonical vertex position. EDM4hep files
    # can store children from one vertex with sub-micron rounding differences
    # in their production positions.
    children_by_parent_set = {}
    for child_index in source_indices:
        child_index = int(child_index)
        parent_set = parents[child_index]
        if parent_set:
            children_by_parent_set.setdefault(parent_set, []).append(child_index)

    assigned_parents = set()
    for parent_set, child_indices in children_by_parent_set.items():
        if assigned_parents.intersection(parent_set):
            raise ValueError(
                f"parents {parent_set} would belong to multiple end vertices"
            )
        assigned_parents.update(parent_set)
        first_parent = parent_set[0]
        parent_endpoint = np.array(
            [
                float(endpoint_x[first_parent]),
                float(endpoint_y[first_parent]),
                float(endpoint_z[first_parent]),
            ]
        )

        # The child whose stored production position is closest to the parent
        # endpoint also supplies the vertex time. Other children can carry
        # small, correlated position/time rounding shifts.
        vertex_distances = [
            np.linalg.norm(
                np.array(
                    [
                        float(vertex_x[child_index]),
                        float(vertex_y[child_index]),
                        float(vertex_z[child_index]),
                    ]
                )
                - parent_endpoint
            )
            for child_index in child_indices
        ]
        reference_child = child_indices[int(np.argmin(vertex_distances))]
        position = (
            float(parent_endpoint[0]),
            float(parent_endpoint[1]),
            float(parent_endpoint[2]),
            float(time[reference_child]) * SPEED_OF_LIGHT_MM_PER_NS,
        )

        # Multiple incoming parents must agree on their common endpoint.
        for parent_index in parent_set[1:]:
            other_parent_endpoint = (
                float(endpoint_x[parent_index]),
                float(endpoint_y[parent_index]),
                float(endpoint_z[parent_index]),
            )
            if not np.allclose(
                other_parent_endpoint, parent_endpoint, rtol=0.0, atol=1.0e-6
            ):
                raise ValueError(
                    f"parents {parent_set} have different endpoint positions"
                )

        vertex = pyhepmc.GenVertex(pyhepmc.FourVector(*position))
        event.add_vertex(vertex)
        for parent_index in parent_set:
            vertex.add_particle_in(particles[parent_index])
        for child_index in child_indices:
            vertex.add_particle_out(particles[child_index])

    # Leave parentless particles at HepMC's implicit root at (0, 0, 0, 0).
    # This extractor does not build explicit source vertices, so reject
    # nonzero source positions/times rather than silently discarding them.
    for particle_index in source_indices:
        particle_index = int(particle_index)
        if parents[particle_index]:
            continue
        position = (
            float(vertex_x[particle_index]),
            float(vertex_y[particle_index]),
            float(vertex_z[particle_index]),
            float(time[particle_index]) * SPEED_OF_LIGHT_MM_PER_NS,
        )
        if not np.allclose(position, 0.0, rtol=0.0, atol=1.0e-9):
            raise ValueError(
                f"parentless particle {particle_index} has nonzero production "
                f"position {position}; explicit source vertices are not "
                "supported by this extractor"
            )

    return event, len(pdg), len(source_indices)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", type=Path, help="npsim EDM4hep ROOT file")
    parser.add_argument("output_file", type=Path, help="output HepMC3 ASCII file")
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="extract only the first N events (useful for a smoke test)",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    input_file = args.input_file.resolve()
    output_file = args.output_file.resolve()

    if not input_file.is_file():
        raise FileNotFoundError(f"input file does not exist: {input_file}")
    if output_file.exists():
        raise FileExistsError(f"refusing to overwrite existing file: {output_file}")
    if not output_file.parent.is_dir():
        raise FileNotFoundError(
            f"output directory does not exist: {output_file.parent}"
        )
    if args.max_events is not None and args.max_events <= 0:
        raise ValueError("--max-events must be a positive integer")

    partial_file = output_file.with_name(output_file.name + ".partial")
    if partial_file.exists():
        raise FileExistsError(f"refusing to overwrite partial file: {partial_file}")

    # Every event is read once; disabling the array cache keeps memory bounded
    # by one event even for large background samples.
    with uproot.open(input_file, array_cache=None) as root_file:
        events = root_file["events"]
        number_of_events = events.num_entries
        if args.max_events is not None:
            number_of_events = min(number_of_events, args.max_events)

        total_particles = 0
        total_source_particles = 0
        # Reserve the partial filename exclusively. Publish the final filename
        # only after the writer closes successfully; retain failures for review.
        with partial_file.open("x"):
            pass
        try:
            with pyhepmc.open(partial_file, "w", format="hepmc3") as writer:
                for entry_number in range(number_of_events):
                    arrays = events.arrays(
                        list(BRANCHES.values()),
                        entry_start=entry_number,
                        entry_stop=entry_number + 1,
                        library="ak",
                    )
                    event, particle_count, source_count = build_hepmc_event(
                        arrays, entry_number
                    )
                    writer.write(event)
                    total_particles += particle_count
                    total_source_particles += source_count

                    print(
                        f"\rExtracted event {entry_number + 1}/{number_of_events}",
                        end="",
                        flush=True,
                    )
            if output_file.exists():
                raise FileExistsError(f"refusing to overwrite: {output_file}")
            partial_file.rename(output_file)
        except Exception:
            print(f"\nIncomplete output retained at {partial_file}", file=sys.stderr)
            raise

    removed_particles = total_particles - total_source_particles
    print(
        f"\nWrote {number_of_events} events and {total_source_particles:,} "
        f"source particles to {output_file}"
    )
    print(
        f"Excluded {removed_particles:,} particles with "
        "BITCreatedInSimulation set"
    )


if __name__ == "__main__":
    main()
