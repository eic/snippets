#!/bin/sh

OPEN=0
PAGES=0

while [ $# -gt 0 ]; do
    case "$1" in
        --open)  OPEN=1;  shift ;;
        --pages) PAGES=1; shift ;;
        *)       break ;;
    esac
done

project="$1"
pipeline="$2"

usage() {
    echo "Usage: $0 [--open] [--pages] <project-id> <pipeline-id>"
    echo ""
    echo "  --open   Serve trace.json locally and open it in Perfetto UI"
    echo "  --pages  Upload trace to eic/perfetto-launcher GitHub Pages and open the launcher URL"
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

# Open trace in Perfetto via a local CORS-enabled HTTP server
open_in_perfetto() {
    _port=$(python3 -c "import socket; s=socket.socket(); s.bind(('',0)); print(s.getsockname()[1]); s.close()")
    _pyserver=$(mktemp /tmp/cors_server_XXXXXX.py)
    # The server runs with no stdin reads; the shell controls shutdown via kill
    cat > "$_pyserver" << 'PYEOF'
import http.server, sys, threading

class CORSHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        super().end_headers()
    def log_message(self, *a):
        pass

port = int(sys.argv[1])
server = http.server.HTTPServer(('127.0.0.1', port), CORSHandler)
server.serve_forever()
PYEOF
    python3 "$_pyserver" "$_port" &
    _server_pid=$!
    # Ensure cleanup even on Ctrl-C or unexpected exit
    trap 'kill "$_server_pid" 2>/dev/null; rm -f "$_pyserver"' EXIT INT TERM
    _launch_url="https://ui.perfetto.dev/#!/?url=http://localhost:${_port}/trace.json"
    echo "Opening: $_launch_url"
    case "$(uname)" in
        Darwin) open "$_launch_url" ;;
        Linux)  xdg-open "$_launch_url" 2>/dev/null || true ;;
        *)      echo "Visit: $_launch_url" ;;
    esac
    echo "Press Enter to stop the local server..."
    read _dummy
    kill "$_server_pid" 2>/dev/null
    rm -f "$_pyserver"
    trap - EXIT INT TERM
}

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

if [ "$OPEN" = "1" ] && [ "$PAGES" = "1" ]; then
    echo "ERROR: --open and --pages are mutually exclusive" >&2
    usage
    exit 1
elif [ "$OPEN" = "1" ]; then
    open_in_perfetto
elif [ "$PAGES" = "1" ]; then
    upload_to_pages "$pipeline"
fi
