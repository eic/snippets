#!/bin/sh

project="$1"
pipeline="$2"

[ -z "$project" -o -z "$pipeline" ] && echo "Usage: $0 <project-id> <pipeline-id>" && exit 1

url="https://eicweb.phy.anl.gov/api/v4/projects"

dir=$(mktemp -d)

print_cmd() {
    echo "$@"
    "$@"
}

fetch_pipeline() {
    project="${1%%:*}"
    pipeline="${1##*:}"
    jobs_url="$url/$project/pipelines/$pipeline/jobs"
    per_page=100
    page=1

    while true; do
        file="$dir/jobs-$project-$pipeline-$page.json"
        # Fetch if the file doesn't exist
        [ -f "$dir/$file" ] || print_cmd curl -LfsS "$jobs_url?include_retried=true&per_page=$per_page&page=$page" -o "$file" || break
        jq -e "length < $per_page" "$file" > /dev/null && break
        page=$((page + 1))
    done

    # Get the bridges, shouldn't be paginated
    bridges_url="$url/$project/pipelines/$pipeline/bridges"

    # Obtain the child pipeline ids
    child_pipelines=$(curl -LfsS "$bridges_url" | jq -r '.[].downstream_pipeline | "\(.project_id):\(.id)"')

    # Fetch child pipelines
    for child_pipeline in $child_pipelines; do
        echo "Fetching pipeline $child_pipeline"
        fetch_pipeline "$child_pipeline"
    done

}

echo "Fetching main pipeline"
fetch_pipeline "$project:$pipeline"

# Finally output as trace.json
jq \
'map(['\
'select(.started_at and .finished_at) | '\
'{name: (.name), cat: "PERF", ph: "B", pid: .pipeline.id, tid: .id, ts: (.started_at | sub("\\.[0-9]+Z$"; "Z") | fromdate * 10e5)},'\
'{name: (.name), cat: "PERF", ph: "E", pid: .pipeline.id, tid: .id, ts: (.finished_at | sub("\\.[0-9]+Z$"; "Z") | fromdate * 10e5)}'\
']) | flatten(1) | .[]' $dir/jobs-*-*.json | jq -s > trace.json

echo "Now upload trace.json to https://ui.perfetto.dev"
