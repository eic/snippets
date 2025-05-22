#!/bin/bash
# This script sets up the environment for EIC simulations using EPIC and EICrecon.
# It installs the necessary software, compiles custom configurations, and generates job submission scripts.
# It also handles the simulation and analysis of events based on user-defined parameters.

# example of usage in BASH(!) shell  - not native .csh shell:
# MOMENTUM=0.5 ./submit.sh

# The workflow includes:
# 1. Setting up the eic-shell
# 2. Compiling EPIC+EICrecon in your directory if not already done.
# 3. Generating a main simulation script = job with ddsim + eicrecon + run root macro.
# 4. Generating a job submission script for HTCondor and submitting it.

#================================================================
#   CONFIGURATION
#================================================================

# simulation parameters
DEFAULT_MOMENTUM=1          # 1 GeV
DEFAULT_PHI=45              # 45 degrees
DEFAULT_THETA=170           # 170 degrees
DEFAULT_PARTICLE=neutron    # "neutron" or "proton"
DEFAULT_NUMBER_OF_EVENTS=10 # 10 events per job
DEFAULT_JOBS=10             # 10 job

# simulation parameters - can be set via command line or environment variables
MOMENTUM=${MOMENTUM:-$DEFAULT_MOMENTUM}
PHI=${PHI:-$DEFAULT_PHI}
THETA=${THETA:-$DEFAULT_THETA}
PARTICLE=${PARTICLE:-$DEFAULT_PARTICLE}
NUMBER_OF_EVENTS=${NUMBER_OF_EVENTS:-$DEFAULT_NUMBER_OF_EVENTS}
JOBS=${JOBS:-$DEFAULT_JOBS}

#================================================================
DETECTOR_CONFIG="epic_backward_hcal_only.xml" # or some other e.g. "epic_backward_hcal_only.xml"
#  make unique dir name for whole configuration
SIM_CONFIG="${PARTICLE}_p${MOMENTUM}gev_phi${PHI}_theta${THETA}_${NUMBER_OF_EVENTS}events"

# output directories
output_dir="/gpfs02/eic/${USER}/output/${SIM_CONFIG}"
my_epic_dir="/gpfs02/eic/${USER}/epic"
my_eicrecon_dir="/gpfs02/eic/${USER}/EICrecon"
current_dir="$(pwd)"

#================================================================
#   UTILITY FUNCTIONS
#================================================================
set -euo pipefail # Exit on error, undefined variables, and pipe failures
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}
error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
    echo "  Context: Function ${FUNCNAME[1]}, Line ${BASH_LINENO[0]}" >&2
}
display_configuration() {
    cat <<EOF
===  SIMULATION CONFIGURATION ===

Simulation Parameters:
  - Particle: ${PARTICLE}
  - Momentum: ${MOMENTUM} GeV
  - Phi: ${PHI}°, Theta: ${THETA}°
  - Events per Job: ${NUMBER_OF_EVENTS}
  - Number of Jobs: ${JOBS}

Output Directory: ${output_dir}
===================================
EOF
}
#=================================================================
# 1. Setting up the eic-shell
# Setup EIC shell environment
# Globals:
#   EICSHELL - Path to eic-shell executable
#   USER - Current username
# Returns:
#   0 on success, 1 on failure
#================================================================
setup_eicshell() {
    export EICSHELL="/eic/u/${USER}/eic/eic-shell"
    if [[ ! -f "$EICSHELL" ]]; then
        log "Installing eic-shell..."
        mkdir -p "$HOME/eic"
        cd "$HOME/eic"
        if ! curl -L https://github.com/eic/eic-shell/raw/main/install.sh | bash; then
            error "Failed to install eic-shell"
            exit 1
        fi
        cd "$current_dir"
    fi
}
#=================================================================

#=================================================================
# 2. Compiling EPIC in your directory if not already done.
#=================================================================

setup_epic() {
    log "Checking EPIC setup..."
    # Check if EPIC is already properly installed
    if [[ -d "$my_epic_dir" ]] && [[ -f "$my_epic_dir/install/bin/thisepic.sh" ]]; then
        log "EPIC already installed at $my_epic_dir"
        return 0
    fi

    log "Setting up EPIC..."

    # Create parent directory if it doesn't exist
    local parent_dir="$(dirname "$my_epic_dir")"
    mkdir -p "$parent_dir"

    # Remove any incomplete installation
    if [[ -d "$my_epic_dir" ]]; then
        log "Removing incomplete EPIC installation..."
        rm -rf "$my_epic_dir"
    fi

    # Clone EPIC
    cd "$parent_dir"
    log "Cloning EPIC repository..."
    if ! git clone https://github.com/eic/epic.git "$(basename "$my_epic_dir")"; then
        error "Failed to clone EPIC repository"
        exit 1
    fi

    cd "$my_epic_dir"
    # Build EPIC

    cat <<EOF | $EICSHELL
    echo "Configuring EPIC build..."

    # Run CMake configure step with ccache launcher support
    echo " Building and installing EPIC (this may take several minutes)..."
    if ! cmake -B build -S . \
        -DCMAKE_INSTALL_PREFIX=install; then
        echo "CMake configuration failed."
        exit 1
    fi
    if ! cmake --build build -j$(nproc) -- install; then
        exit 1
    fi
EOF
    # Verify installation
    if [[ ! -f "$my_epic_dir/install/bin/thisepic.sh" ]]; then
        error "EPIC installation incomplete - missing thisepic.sh"
        exit 1
    fi

    cd "$current_dir"
    log "EPIC setup completed successfully"

}

#================================================================
#   3. Compiling EICrecon in your directory if not already done.
#================================================================

setup_eicrecon() {
    log "Checking EICrecon setup..."
    # First, ensure EPIC is properly installed
    if [[ ! -f "$my_epic_dir/install/bin/thisepic.sh" ]]; then
        log "Expected EPIC installation at: $my_epic_dir/install/bin/thisepic.sh"
        error "EPIC installation not found. Please set up EPIC first."
        exit 1
    fi
    # Check if EICrecon is already properly installed
    if [[ -d "$my_eicrecon_dir" ]] && [[ -f "$my_eicrecon_dir/install/bin/eicrecon-this.sh" ]]; then
        log "EICrecon already installed at $my_eicrecon_dir"
        return 0
    fi

    log "Setting up EICrecon..."

    # Create parent directory if it doesn't exist
    local parent_dir="$(dirname "$my_eicrecon_dir")"
    mkdir -p "$parent_dir"

    # Remove any incomplete installation
    if [[ -d "$my_eicrecon_dir" ]]; then
        log "Removing incomplete EICrecon installation..."
        rm -rf "$my_eicrecon_dir"
    fi

    # Clone EICrecon
    cd "$parent_dir"
    log "Cloning EICrecon repository..."
    if ! git clone https://github.com/eic/EICrecon.git "$(basename "$my_eicrecon_dir")"; then
        error "Failed to clone EICrecon repository"
        exit 1
    fi

    # Build EICrecon
    cd "$my_eicrecon_dir"

    log "Configuring EICrecon build..."
    cat <<EOF | $EICSHELL
    source "$my_epic_dir/install/bin/thisepic.sh"
    if ! cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install; then
        exit 1
    fi
    echo "Building and installing EICrecon (this may take several minutes)..."
    if ! cmake --build build -j$(nproc) -- install; then
        exit 1
    fi
    echo "Installing EICrecon..."
    if ! cmake --install build; then
        exit 1
    fi
EOF
    # Verify installation
    if [[ ! -f "$my_eicrecon_dir/install/bin/eicrecon-this.sh" ]]; then
        error "EICrecon installation incomplete - missing eicrecon-this.sh"
        exit 1
    fi

    cd "$current_dir"
    log "EICrecon setup completed successfully"
}
#=================================================================

#================================================================
#   MAIN SCRIPT GENERATION
#================================================================

generate_main_script() {
    log "Generating main simulation script..."

    cat >my_generated_script.sh <<'EOF'
#!/bin/bash
cd OUTPUT_DIR_PLACEHOLDER
# Parse arguments
if [[ $# -ne 7 ]]; then
    echo "Usage: $0 <Cluster> <Process> <Momentum> <Phi> <Theta> <Particle> <NumberOfEvents>"
    echo "current values: $CLUSTER $PROCESS $GUN_MOMENTUM $GUN_PHI $GUN_THETA $PARTICLE $NUMBER_OF_EVENTS"
    exit 1
fi
export CLUSTER="$1"
export PROCESS="$2"
export GUN_MOMENTUM="$3"
export GUN_PHI="$4"
export GUN_THETA="$5"
export PARTICLE="$6"
export NUMBER_OF_EVENTS="$7"

# Set up file names
export FILENAME="${PARTICLE}_${NUMBER_OF_EVENTS}events_p${GUN_MOMENTUM}gev_phi${GUN_PHI}_theta${GUN_THETA}_job${CLUSTER}_${PROCESS}"
export DDSIM_FILE="sim_${FILENAME}.edm4hep.root"
export EICRECON_FILE="eicrecon_${FILENAME}.edm4eic.root"

# Calculate angle ranges (small ranges for single-angle shooting)
export GUN_THETA_MIN=$(echo "$GUN_THETA - 0.0001" | bc -l)
export GUN_THETA_MAX=$(echo "$GUN_THETA + 0.0001" | bc -l)
export GUN_PHI_MIN=$(echo "$GUN_PHI - 0.0001" | bc -l)
export GUN_PHI_MAX=$(echo "$GUN_PHI + 0.0001" | bc -l)
export GUN_MOMENTUM_MIN=$(echo "$GUN_MOMENTUM - 0.00001" | bc -l)
export GUN_MOMENTUM_MAX=$(echo "$GUN_MOMENTUM + 0.00001" | bc -l)

cat << 'EOFINNER' | EICSHELL_PLACEHOLDER
    # Source environment
    if [[ -f "MY_EPIC_DIR_PLACEHOLDER/install/bin/thisepic.sh" ]]; then
        source "MY_EPIC_DIR_PLACEHOLDER/install/bin/thisepic.sh" epic
    else
        error "EPIC installation not found"
        exit 1
    fi

    if ! ddsim \
            --compactFile "$DETECTOR_PATH/DETECTOR_CONFIG_PLACEHOLDER" \
            --numberOfEvents "$NUMBER_OF_EVENTS" \
            --random.seed "$(date +%N)" \
            --enableGun \
            --gun.particle "$PARTICLE" \
            --gun.thetaMin "${GUN_THETA_MIN}*degree" \
            --gun.thetaMax "${GUN_THETA_MAX}*degree" \
            --gun.phiMin "${GUN_PHI_MIN}*degree" \
            --gun.phiMax "${GUN_PHI_MAX}*degree" \
            --gun.distribution uniform \
            --gun.momentumMin "${GUN_MOMENTUM_MIN}*GeV" \
            --gun.momentumMax "${GUN_MOMENTUM_MAX}*GeV" \
            --outputFile "$DDSIM_FILE"; then
            echo "DDSIM simulation failed"
            exit 1
    fi

    # Source EICrecon
    if [[ -f "MY_EICRECON_DIR_PLACEHOLDER/install/bin/eicrecon-this.sh" ]]; then
        source "MY_EICRECON_DIR_PLACEHOLDER/install/bin/eicrecon-this.sh" epic
    else
        echo "EICrecon installation not found"
        exit 1
    fi
    # Run EICrecon if needed
    echo "Running EICrecon..."
    if ! eicrecon "$DDSIM_FILE" \
            -Ppodio:output_file="$EICRECON_FILE"
            -Ppodio:output_collections="MCParticles"; then
        echo "EICrecon failed"
        exit 1
    fi
    # Run analysis
    echo "Running ROOT analysis..."
    analysis_script="CURRENT_DIR_PLACEHOLDER/example_macro.C"
    if [[ ! -f "$analysis_script" ]]; then
        echo "Analysis script not found: $analysis_script"
        exit 1
    fi
    output_file="ana_${FILENAME}.root"
    if ! root -l -b -q "${analysis_script}(\\\"${EICRECON_FILE}\\\", \\\"${output_file}\\\")"; then
        echo "ROOT analysis failed"
        exit 1
    fi
    echo "Job completed successfully"
EOFINNER
EOF
    # Replace placeholders with actual values
    sed -i "s|OUTPUT_DIR_PLACEHOLDER|$output_dir|g" my_generated_script.sh
    sed -i "s|MY_EPIC_DIR_PLACEHOLDER|$my_epic_dir|g" my_generated_script.sh
    sed -i "s|MY_EICRECON_DIR_PLACEHOLDER|$my_eicrecon_dir|g" my_generated_script.sh
    sed -i "s|CURRENT_DIR_PLACEHOLDER|$current_dir|g" my_generated_script.sh
    sed -i "s|EICSHELL_PLACEHOLDER|$EICSHELL|g" my_generated_script.sh
    sed -i "s|DETECTOR_CONFIG_PLACEHOLDER|$DETECTOR_CONFIG|g" my_generated_script.sh

    chmod +x my_generated_script.sh

    log "Main simulation script generated successfully"
}

#================================================================
#   JOB SUBMISSION SCRIPT
#================================================================
generate_job_script() {
    local temp_job="$current_dir/generated.job"
    # Remove existing job file if it exists
    if [[ -f "$temp_job" ]]; then
        rm -f "$temp_job"
    fi
    cat >"$temp_job" <<EOF
Universe                = vanilla
GetEnv                  = False
Requirements            = (CPU_Speed >= 1)
Rank                    = CPU_Speed
Initialdir              = $current_dir
Arguments               = \$(Cluster) \$(Process) $MOMENTUM $PHI $THETA $PARTICLE $NUMBER_OF_EVENTS
Executable              = my_generated_script.sh
Error                   = $output_dir/log/error\$(Cluster)_\$(Process).err
Output                  = $output_dir/log/out\$(Cluster)_\$(Process).out
Log                     = $output_dir/log/log\$(Cluster)_\$(Process).log
Queue $JOBS
EOF
    if [[ ! -f "$temp_job" ]]; then
        return 1
    fi
    # Only echo the filename - no log messages
    echo "$temp_job"
}

#================================================================
#   MAIN EXECUTION
#================================================================
main() {
    log "=== Starting EPIC Simulation ==="
    #================================================================
    display_configuration

    # Setup phase
    log "Setting up eic-shell..."
    setup_eicshell
    setup_epic
    setup_eicrecon

    # Create output directories
    mkdir -p "$output_dir/log"
    rm -f $output_dir/log/*.*
    log "Output directory created: $output_dir"

    generate_main_script

    log "Generating HTCondor job submission script..."
    local job_file=$(generate_job_script)

    if [[ -z "$job_file" ]] || [[ ! -f "$job_file" ]]; then
        error "Failed to create job submission file"
        exit 1
    fi
    log "Job submission file created: $job_file"

    log "Submitting job..."
    if condor_submit "$job_file"; then
        log "Job submitted successfully! Output files will be saved in $output_dir"
    else
        error "Job submission failed"
        exit 1
    fi
    # wait for the jobs to finish
    #log "Waiting for jobs to finish..."
    # ./condor_control.sh
    log "=== EPIC Simulation Completed ==="
    #  now one can merge the output files into one
    # hadd -f -j -k "${output_dir}/merged_ana.root" "${output_dir}/ana*.root"
}
# Run main function
main "$@"
