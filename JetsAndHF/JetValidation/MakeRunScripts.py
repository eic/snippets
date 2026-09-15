# =============================================================================
#! @file    MakeRunScripts.py
#! @authors Derek Anderson (derek.murphy.anderson@protonmail.com)
# -----------------------------------------------------------------------------
#! @brief Helper script to generate various scripts to hand-off
#!   to eic-shell. This is to get around arguments not being
#!   able to be passed to scripts ran in eic-shell.  
#!
#! @usage
#!     MakeRunScripts.py -d <my rucio DID> \
#!                       -c <descriptive tag for campaign> \
#!                       -a <descriptive tag for dataset> \
#!                       -s <number of lists to split files into> \
#!                       -l <padding length for indices>
# =============================================================================

from typing import List
import argparse as ap
import os

class ScriptMaker:
    def __init__(self, did: str, camp: str, data: str, splits: int, length: int):
        self.did    = did
        self.camp   = camp
        self.data   = data
        self.splits = splits
        self.length = length
        self.suffix = f"epic{self.camp}_{self.data}"

    def _make_split_names(self, pre: str, ext: str) -> List[str]:
        return [f"{pre}.split_{split:0{self.length}}.{self.suffix}.{ext}" for split in range(self.splits)]

    def make_rucio_script(self) -> None:
        script_name = f"run/scripts/do_rucio.{self.suffix}.sh"
        list_name   = f"run/lists/files.{self.suffix}.list"
        split_pref  = f"run/lists/files.split_"
        list_comm   = f"rucio replica list file --protocols root --pfns --rses isopenaccess {self.did} > {list_name}"
        split_comm  = f"split --suffix-length={self.length} --additional-suffix=.{self.suffix}.list -d -n l/{self.splits} {list_name} {split_pref}"

        # generate script
        os.makedirs(os.path.dirname(script_name), exist_ok=True)
        os.makedirs(os.path.dirname(list_name), exist_ok=True)
        with open(script_name, 'w') as script:
            script.write("#!/bin/bash\n\n")
            script.write(f"{list_comm}\n")
            script.write(f"{split_comm}\n")
        os.chmod(script_name, 0o777)

    def make_hist_scripts(self) -> None:
        drivers   = self._make_split_names("run/scripts/do_hists", "sh")
        in_lists  = self._make_split_names("run/lists/files", "list")
        out_files = self._make_split_names("run/hists/hists", "root")
        for split, driver in enumerate(drivers):

            # generate command
            hist_comm = f'root -b -q "MakeJetValidationHists.C(\\\"{out_files[split]}\\\", \\\"{in_lists[split]}\\\", 1)"'

            # generate script
            os.makedirs(os.path.dirname(out_files[split]), exist_ok=True)
            with open(driver, 'w') as script:
                script.write("#!/bin/bash\n\n")
                script.write(hist_comm)
            os.chmod(driver, 0o777)

    def make_plot_scripts(self) -> None:
        hist_files  = " ".join(self._make_split_names("run/hists/hists", "root"))
        merge_file  = f"run/hists/hist.merged.{self.suffix}.root"
        merge_comm  = f"hadd {merge_file} {hist_files}"
        plot_comm   = f'root -b -q "MakeJetValidationPlots.C(\\\"run/plots\\\", \\\"{self.suffix}\\\", \\\"{merge_file}\\\")'
        script_name = f"run/scripts/do_plots.{self.suffix}.sh"

        # generate script
        os.makedirs("run/plots", exist_ok=True)
        with open(script_name, 'w') as script:
            script.write("#!/bin/bash\n\n")
            script.write(f"{merge_comm}\n")
            script.write(f"{plot_comm}\n")
        os.chmod(script_name, 0o777)


if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument("-d", "--did")
    parser.add_argument("-c", "--camp")
    parser.add_argument("-a", "--data")
    parser.add_argument("-s", "--splits")
    parser.add_argument("-l", "--length")
    args = parser.parse_args()

    maker = ScriptMaker(args.did, args.camp, args.data, int(args.splits), int(args.length))
    maker.make_rucio_script()
    maker.make_hist_scripts()
    maker.make_plot_scripts()
