#!/usr/bin/env python3
"""Deterministic physics/IO checks; no network, ROOT bindings, or pytest needed.

Run: python test_svt_rate.py. Generated outputs are retained in ./tmp by default;
set SVT_TEST_OUTPUT to choose another validation output directory.
"""

import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import awkward as ak
import numpy as np

import plot_hit_rate_per_source as rate
from svt_sources import RELATIONS, assign_sources, read_branch, resolve_status


class Branch:
    def __init__(self, data):
        self.data = data

    def array(self, library, entry_start, entry_stop):
        return self.data[entry_start:entry_stop]


class RootFile:
    def __init__(self, tree):
        self.tree = tree

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def __getitem__(self, name):
        assert name == "events"
        return self.tree


class Tree(dict):
    num_entries = 2


def fixture():
    """Two equal-duration frames, including one empty frame in normalization."""
    tree = Tree()

    def add(name, fields):
        tree[name] = Branch(ak.zip({name + "." + field: ak.Array([values, []])
                                   for field, values in fields.items()}))

    add("MCParticles", {"generatorStatus": [1, 0, 2001, 3001, 9999],
                        "parents_begin": [0, 0, 1, 1, 1], "parents_end": [0, 1, 1, 1, 1]})
    add("_MCParticles_parents", {"index": [0]})
    for config in RELATIONS:
        if config["rec"] == "SiBarrelVertexRecHits":
            coordinates = {"position.x": [35.] * 8, "position.y": [0.] * 8, "position.z": [0.] * 8}
            rec_raw = [0, 1, 2, 3, 4, 5, 6, -1]
            raw, sim, weights, particles = [0, 1, 2, 2, 3, 4], [0, 1, 2, 3, -1, 4], [.5, 1., .2, .8, 1., 1.], [0, 1, 0, 2, 4]
        else:
            endcap = config["rec"] == "SiEndcapTrackerRecHits"
            coordinates = {"position.x": [10. if endcap else 270.], "position.y": [0.], "position.z": [250. if endcap else 0.]}
            rec_raw, raw, sim, weights, particles = [0], [0], [0], [1.], [2 if endcap else 3]
        add(config["rec"], coordinates)
        add(config["rec_raw"], {"index": rec_raw})
        add(config["association"], {"weight": weights})
        add(config["association_raw"], {"index": raw})
        add(config["association_sim"], {"index": sim})
        add(config["sim_particle"], {"index": particles})
    return tree


def args(chunk_size=1):
    return SimpleNamespace(tile=False, entry_start=0, entry_stop=None,
                           chunk_size=chunk_size, frame_us=2.)


class RateTests(unittest.TestCase):
    def test_status_tracing_and_invalid_ancestry(self):
        self.assertEqual(resolve_status(1, [1, 0], [0, 0], [0, 1], [0]), (1, True, False))
        self.assertEqual(resolve_status(0, [0], [0], [1], [0]), (0, True, True))
        self.assertEqual(resolve_status(-1, [1], [0], [0], []), (None, False, False))
        self.assertEqual(resolve_status(0, [0], [0], [2], [0]), (0, True, True))

    def test_weight_choice_secondary_unmatched_invalid_and_other(self):
        tree = fixture()
        config = RELATIONS[1]
        data = {key: read_branch(tree, value, 0, 1)[0] for key, value in config.items()}
        labels, diagnostic = assign_sources(data["rec_raw"], data["association"]["weight"],
                                            data["association_raw"], data["association_sim"],
                                            data["sim_particle"], read_branch(tree, "MCParticles", 0, 1)[0],
                                            read_branch(tree, "_MCParticles_parents", 0, 1)[0])
        self.assertEqual(labels.tolist(), ["DIS", "DIS", "SR", "invalid_relation", "other", "unmatched", "unmatched", "invalid_relation"])
        self.assertEqual(diagnostic["secondary_traced_count"], 1)
        self.assertEqual(diagnostic["ambiguous_source_count"], 1)
        self.assertEqual(diagnostic["invalid_relation_count"], 2)
        bad_weights = ak.Array([float("nan")] * 6)
        labels, diagnostic = assign_sources(data["rec_raw"], bad_weights, data["association_raw"], data["association_sim"], data["sim_particle"], read_branch(tree, "MCParticles", 0, 1)[0], read_branch(tree, "_MCParticles_parents", 0, 1)[0])
        self.assertEqual(diagnostic["invalid_weight_count"], 5)

    def test_chunks_conservation_zero_frames_and_fractional_rate(self):
        for tile in (False, True):
            results = []
            for size in (1, 2):
                options = args(size)
                options.tile = tile
                with patch.object(rate.uproot, "open", side_effect=lambda *a, **k: RootFile(fixture())):
                    maps, inputs, diagnostics, events, elapsed = rate.process_files(["synthetic.root"], options)
                self.assertEqual(events, 2)
                self.assertEqual(inputs[0]["processed_events"], 2)
                self.assertEqual(maps["L0"]["counts"]["ALL"].sum(), 8)
                np.testing.assert_array_equal(maps["L3"]["counts"]["other background"], maps["L3"]["counts"]["Bremstrahlung"])
                results.append(maps)
            for name in results[0]:
                for source in results[0][name]["counts"]:
                    np.testing.assert_array_equal(results[0][name]["counts"][source], results[1][name]["counts"][source])
        # 3 counts over 3000 frames of 2 us => 0.5 kHz; never truncate to zero.
        rows = rate.summarize(results[0], 1000 / (2 * 3000))
        dis = next(row for row in rows if row["lname"] == "L0" and row["source"] == "DIS")
        self.assertAlmostEqual(dis["rsu_or_tile_max_rate_khz"], 1 / 3)

    def test_map_edge_coverage_and_outside_counter(self):
        maps = rate.make_surface_maps()
        hits = ak.Array([{"position.x": 35., "position.y": 0., "position.z": 195.},
                         {"position.x": 35., "position.y": 0., "position.z": 205.}])
        rate.fill_event(maps, "SiBarrelVertexRecHits", hits, np.array(["DIS", "DIS"]))
        self.assertEqual(maps["L0"]["counts"]["ALL"].sum(), 1)
        self.assertEqual(maps["L0"]["counters"]["ALL"]["outside_map_hits"], 1)
        self.assertEqual(maps["L0"]["xedges"][-1], 200.)

    def test_missing_fields_and_relation_length(self):
        with self.assertRaises(ValueError):
            read_branch(fixture(), "MCParticles", 0, 1, ["missing"])
        with self.assertRaises(ValueError):
            rate.fill_event(rate.make_surface_maps(), "SiBarrelVertexRecHits", ak.Array([]), np.array(["DIS"]))

    def test_cli_outputs_and_input_list(self):
        parent = Path(os.environ.get("SVT_TEST_OUTPUT", str(rate.HERE / "tmp")))
        parent.mkdir(parents=True, exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix="svt_rate_validation_", dir=str(parent)))
        listing = output / "inputs.txt"
        listing.write_text("# test fixture\n\nsynthetic.root:events\n")
        inputs = rate.build_inputs(SimpleNamespace(file_list=listing, input_files=None, max_files=None))
        self.assertEqual(inputs, [str(output / "synthetic.root")])
        command = ["plot_hit_rate_per_source.py", "--file-list", str(listing), "--output-dir", str(output), "--label", "smoke", "--chunk-size", "1"]
        with patch.object(sys, "argv", command), patch.object(rate.uproot, "open", side_effect=lambda *a, **k: RootFile(fixture())):
            rate.main()
        result = output / "smoke"
        record = json.loads((result / "smoke_run.json").read_text())
        self.assertEqual(record["processed_events"], 2)
        self.assertEqual(record["rate_scale_khz_per_count"], 250.)
        with (result / "smoke_summary.csv").open() as stream:
            rows = list(csv.DictReader(line for line in stream if not line.startswith("#")))
        all_l0 = next(row for row in rows if row["lname"] == "L0" and row["source"] == "ALL")
        self.assertEqual(float(all_l0["scaled_total_rate_khz"]), 2000.)
        with np.load(result / "smoke_maps.npz") as archive:
            index = next(row for row in record["maps"] if row["lname"] == "L0" and row["source"] == "ALL")
            self.assertEqual(archive[index["key"]].sum(), 8)
        for filename in record["pdfs"]:
            self.assertTrue((result / filename).read_bytes().startswith(b"%PDF"))
        with patch.object(sys, "argv", command), self.assertRaises(FileExistsError):
            rate.main()
        print(f"Validation outputs retained at {output}")


if __name__ == "__main__":
    unittest.main()
