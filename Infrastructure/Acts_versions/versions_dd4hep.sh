git clone https://github.com/AIDASoft/DD4hep.git
git clone https://github.com/eic/containers.git
git -C DD4hep for-each-ref --format="%(refname:short)|%(creatordate:format:%b %d %Y %H:%M:%S)" "refs/tags/*" | tee tags
git -C containers for-each-ref --format="%(refname:short)|%(creatordate:format:%b %d %Y %H:%M:%S)" "refs/tags/*" | grep 'stable' | grep -v '^v23' | grep -v '^v24.02' | tee eic_container_tags
cat eic_container_tags | sed -e 's#^v##' | awk -F'|' '{print $1}' | xargs -I{} singularity exec /cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:{} sh -c "spack find --json dd4hep | jq -r '.[].version'" | sed -e 's/\./\.\\\+/g' | xargs -I {} sh -c "git -C DD4hep tag | grep -E '{}$'" | xargs -I{} git -C DD4hep for-each-ref --format="%(refname:short)|%(creatordate:format:%b %d %Y %H:%M:%S)" "refs/tags/{}" | tee result
