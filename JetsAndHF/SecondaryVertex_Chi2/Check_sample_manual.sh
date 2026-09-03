#!/bin/bash
# Script to check the simulated sample from July campaign
# Shyam Kumar; University of Houston; shyam.kumar@cern.ch 
BLUE='\033[1;34m'
RED='\033[1;31m'
RESET='\033[0m'

# eAg DIS sample
eAgDISSample=$(rucio did list --short "epic:/RECO/26.07.*/epic_craterlake/DIS/BeAGLE*/eAg/9x115/*")

if [[ -n "$eAgDISSample" ]]; then
    echo -e "${BLUE} eAg (9x115) DIS Sample found${RESET}"
    echo -e "${eAuDISSample}\n"

else
    echo -e "${RED}eAg (9x115) DIS Sample not found${RESET}"
fi
# eAu DIS sample
eAuDISSample=$(rucio did list --short "epic:/RECO/26.07.*/epic_craterlake/DIS/BeAGLE*/eAu/9x100/*")

if [[ -n "$eAuDISSample" ]]; then
    echo -e "${BLUE}eAu (9x100) DIS Sample found${RESET}"
    echo -e "${eAuDISSample}\n"
else
    echo -e "${RED}eAu (9x100) DIS Sample not found${RESET}"
fi

# eAu D0, Lc, Ds sample
enriched_sample=("D0" "Lc" "Ds")

for ((i = 0; i < ${#enriched_sample[@]}; i++)); do
    particle="${enriched_sample[$i]}"

    eAuSample=$(rucio did list --short \
        "epic:/RECO/26.07.*/epic_craterlake/SIDIS/${particle}_ABCONV/*/eAu/9x100/*")

    if [[ -n "$eAuSample" ]]; then
        echo -e "${BLUE}eAu (9x100) ${particle} sample found ${RESET}"
        echo -e "${eAuSample}\n"
    else
        echo -e "${RED}eAu (9x100) ${particle} sample not found ${RESET}"
    fi
done
# ep DIS sample (9x130)
epDISSample=$(rucio did list --short "epic:/RECO/26.07.*/epic_craterlake/DIS/pythia8*/NC/noRad/ep/9x130/*")

if [[ -n "$epDISSample" ]]; then
    echo -e "${BLUE}ep (9x130) DIS Sample found ${RESET}"
    echo -e "${epDISSample}\n"
else
    echo -e "${RED}ep (9x130) DIS Sample not found ${RESET}"
fi

# ep D0 sample (9x130)
epD0Sample=$(rucio did list --short "epic:/RECO/26.07.*/epic_craterlake/SIDIS/D0_ABCONV/*/*/ep/9x130/*")

if [[ -n "$epD0Sample" ]]; then
    echo -e "${BLUE}ep (9x130) D0 Sample found ${RESET}"
    echo -e "${epD0Sample}\n"

else
    echo -e "${RED}ep (9x130) D0 Sample not found ${RESET}"
fi

# ep DIS sample (9x275)
epDISSample=$(rucio did list --short "epic:/RECO/26.07.*/epic_craterlake/DIS/pythia8*/NC/noRad/ep/9x275/*")

if [[ -n "$epDISSample" ]]; then
    echo -e "${BLUE}ep (9x275) DIS Sample found ${RESET}"
    echo -e "${epDISSample}\n"
else
    echo -e "${RED}ep (9x275) DIS Sample not found ${RESET}"
fi

for ((i = 0; i < ${#enriched_sample[@]}; i++)); do
    particle="${enriched_sample[$i]}"

    epSample=$(rucio did list --short \
        "epic:/RECO/26.07.*/epic_craterlake/SIDIS/${particle}_ABCONV/*/*/ep/9x275/*")
    if [[ -n "$epSample" ]]; then
        echo -e "${BLUE}ep (9x275) ${particle} sample found ${RESET}"
        echo -e "${epSample}\n"
    else
        echo -e "${RED}ep (9x275) ${particle} sample not found ${RESET}"
    fi
done

