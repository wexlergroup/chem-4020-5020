#!/bin/bash
# Builds every studio in development/ as one self-contained HTML file and
# copies it to book/computational-studios/, which the Pages deploy publishes.
set -euo pipefail
cd "$(dirname "$0")"
for app in */; do
  app=${app%/}
  [ -f "$app/package.json" ] || continue
  (cd "$app" && npm ci && npx vite build)
  cp "$app/dist/index.html" "../book/computational-studios/$app.html"
done
