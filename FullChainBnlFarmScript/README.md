# Instructions to use for BNL HTCondor farm:

1. Login to eic node
2. Run this command with different parameters (it will use default if not specified)
```bash
MOMENTUM=3 ./run_submit.sh
```
3. You could modify DEFAULT parameters inside the script.
4. Example in loop:
``` bash
momenta=(1 2 5)
for mom in "${momenta[@]}"; do
    export MOMENTUM=$mom
    ./run.sh
done
```
5. This will set everything for you, including `eic-shell`, `epic` repository and `eicrecon`.
And run in batch mode full chain:
```bash
ddsim SIMFILE
eicrecon RECOFILE
root -l -n yourMacro.C (RECOFILE, ANAFILE)
```
6. There is a possibility to wait for jobs using another script `condor_control.sh` which monitors your batch jobs and continues when everything is finished( e.g. you could merge output root analized files which contain your own trees/histograms)

7. The script is a bit long (~400 lines) yet it does fully automated procedure in single file. 
