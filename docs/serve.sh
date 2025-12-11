#!/bin/bash
# Simple HTTP server to view documentation locally
# Usage: ./serve.sh [port]

PORT=${1:-8000}
DOCS_DIR="$(dirname "$0")/build"

if [ ! -d "$DOCS_DIR" ]; then
    echo "Error: Documentation not built. Run: julia --project=docs docs/make.jl"
    exit 1
fi

echo "Starting documentation server..."
echo "Open your browser to: http://localhost:$PORT"
echo "Press Ctrl+C to stop the server"
echo ""

cd "$DOCS_DIR"
python3 -m http.server "$PORT"
