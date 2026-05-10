#!/bin/sh

PAGES=0

while [ $# -gt 0 ]; do
    case "$1" in
        --pages) PAGES=1; shift ;;
        *)       break ;;
    esac
done

project="$1"
pipeline="$2"

usage() {
    echo "Usage: $0 [--pages] <project-id> <pipeline-id>"
    echo ""
    echo "  --pages  Upload trace to eic/perfetto-launcher GitHub Pages and print the launcher URL"
    echo "           Requires GITHUB_PAGES_TOKEN env var with contents:write on eic/perfetto-launcher"
}

[ -z "$project" -o -z "$pipeline" ] && usage && exit 1

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

    # Obtain the child pipeline ids; skip bridges where the child pipeline
    # has not yet been created (downstream_pipeline is null)
    child_pipelines=$(curl -LfsS "$bridges_url" | jq -r '.[] | select(.downstream_pipeline != null) | .downstream_pipeline | "\(.project_id):\(.id)"')

    # Fetch child pipelines
    for child_pipeline in $child_pipelines; do
        echo "Fetching pipeline $child_pipeline"
        fetch_pipeline "$child_pipeline"
    done

}

echo "Fetching main pipeline"
fetch_pipeline "$project:$pipeline"

# Finally output as trace.json; skip jobs with null pipeline id (bridge
# triggered but child pipeline not yet created when waterfall ran)
jq \
'map(['\
'select(.started_at and .finished_at and .pipeline.id != null) | '\
'{name: (.name), cat: "PERF", ph: "B", pid: .pipeline.id, tid: .id, ts: (.started_at | sub("\\.[0-9]+Z$"; "Z") | fromdate * 10e5)},'\
'{name: (.name), cat: "PERF", ph: "E", pid: .pipeline.id, tid: .id, ts: (.finished_at | sub("\\.[0-9]+Z$"; "Z") | fromdate * 10e5)}'\
']) | flatten(1) | .[]' $dir/jobs-*-*.json | jq -s > trace.json

echo "Now upload trace.json to https://ui.perfetto.dev"

# Upload trace to eic/perfetto-launcher GitHub Pages and print the launcher URL
upload_to_pages() {
    _pipeline="$1"
    _dest="traces/trace-${_pipeline}.json"
    _launcher_url="https://eic.github.io/perfetto-launcher/?trace=${_dest}"

    if [ -z "$GITHUB_PAGES_TOKEN" ]; then
        echo "ERROR: GITHUB_PAGES_TOKEN is not set; cannot upload to GitHub Pages" >&2
        return 1
    fi

    echo "Uploading trace to eic/perfetto-launcher as ${_dest}..."

    # Fetch existing file sha if it exists (needed for updates via GitHub Contents API)
    _sha=$(curl -sS \
        -H "Authorization: token ${GITHUB_PAGES_TOKEN}" \
        "https://api.github.com/repos/eic/perfetto-launcher/contents/${_dest}" \
        | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('sha',''))" 2>/dev/null || true)

    # Build the JSON body with Python to avoid shell "Argument list too long" on large traces
    _body_file=$(mktemp /tmp/gh_upload_XXXXXX.json)
    python3 - << PYEOF
import json, base64
with open('trace.json', 'rb') as f:
    content = base64.b64encode(f.read()).decode()
sha = '${_sha}'
body = {
    'message': ('Update' if sha else 'Add') + ' trace for pipeline ${_pipeline}',
    'content': content,
}
if sha:
    body['sha'] = sha
with open('${_body_file}', 'w') as out:
    json.dump(body, out)
PYEOF

    _response=$(curl -sS -X PUT \
        -H "Authorization: token ${GITHUB_PAGES_TOKEN}" \
        -H "Content-Type: application/json" \
        "https://api.github.com/repos/eic/perfetto-launcher/contents/${_dest}" \
        -d "@${_body_file}" \
        -w "\n%{http_code}")
    rm -f "$_body_file"
    _http_code=$(echo "$_response" | tail -1)
    if [ "$_http_code" != "201" ] && [ "$_http_code" != "200" ]; then
        echo "ERROR: GitHub API returned HTTP ${_http_code}" >&2
        echo "$_response" | head -n -1 >&2
        return 1
    fi

    echo "Trace uploaded successfully."
    echo "Launcher URL: ${_launcher_url}"
    echo "(traces/index.json is regenerated automatically by the deploy workflow)"

    if [ -t 1 ]; then
        case "$(uname)" in
            Darwin) open "$_launcher_url" ;;
            Linux)  xdg-open "$_launcher_url" 2>/dev/null || true ;;
        esac
    fi
}

if [ "$PAGES" = "1" ]; then
    upload_to_pages "$pipeline"
fi
