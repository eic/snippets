# Generate material map for ePIC

See run_matmap_py.sh for options and comments. Scripts depend on ACTS (version 30+) python bindings. Tested to work with eic-shell v25.01.0 (ACTS v34.1).
* NOTE: the script experienced a lot of FPE errors with ACTS34. Temporarily disabled the monitor by:
`export ACTS_SEQUENCER_DISABLE_FPEMON=0`

For more info, please check slides at https://indico.bnl.gov/event/22390/

To run EICrecon with a local material map, use flag -Pacts:MaterialMap=XXX 


