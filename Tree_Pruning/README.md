# Tree Pruning Script

This is (appropriately) a trimmed down version of a short script I made for the analysis tutorials. The script takes in an EICrecon output file and trims it down to a smaller subset of branches. This might be useful for a few reasons -

1. Easier for new users to digest a simpler file (avoid getting overwhelmed with the number of branches there by default)
2. Cut down on the size of files if you need to copy them (but don't need all of the information contained within them.

There's some basic checking in terms of whether the file exists (and if it's a root file) at the start. Note that currently, it does NOT check if branches exist before trying to grab them. Maybe one for the future.

## Execution

Execute via -

root TreePrune.C

You will be prompted for an input file. You could specify this at execution too -

root 'TreePrune.C("InputFilePath")'

## Disclaimer

I tested this by pruning files utilised in the analysis tutorial and checking the output from running the exercise scripts on the pruned and unpruned versions. The output was identical. It should probably be tested a bit further if you want to rely on it for something a bit beyond this though.