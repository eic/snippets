git clone https://github.com/acts-project/acts
git clone https://github.com/eic/containers.git
git -C acts for-each-ref --format="%(refname:short)|%(creatordate:format:%b %d %Y %H:%M:%S)" "refs/tags/*" | tee tags
git -C containers for-each-ref --format="%(refname:short)|%(creatordate:format:%b %d %Y %H:%M:%S)" "refs/tags/*" | grep 'stable' | grep -v '^v23' | grep -v '^v24.02' | tee eic_container_tags
cat eic_container_tags | sed -e 's#^v##' | awk -F'|' '{print $1}' | xargs -I{} singularity exec /cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:{} sh -c "spack find --json acts | jq -r '.[].version'" | xargs -I{} git -C acts for-each-ref --format="%(refname:short)|%(creatordate:format:%b %d %Y %H:%M:%S)" "refs/tags/v{}" | tee result