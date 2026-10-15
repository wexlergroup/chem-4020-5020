#!/bin/bash
# Installs the tools needed to build the Jupyter Book in Claude Code cloud
# sessions. Mirrors .github/workflows/deploy.yml: the pinned Python packages
# in book/requirements.txt plus the MyST CLI from npm.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "$CLAUDE_PROJECT_DIR"

VENV="$CLAUDE_PROJECT_DIR/.venv"
if [ ! -x "$VENV/bin/python" ]; then
  python3 -m venv "$VENV"
fi

"$VENV/bin/python" -m pip install --quiet --upgrade pip
"$VENV/bin/python" -m pip install --quiet -r book/requirements.txt
"$VENV/bin/python" -m ipykernel install --sys-prefix --name python3 >/dev/null

# Cloud sessions run as root, and Jupyter Server refuses to start as root by
# default, which makes `myst build --execute` fail to run any code cells.
mkdir -p "$VENV/etc/jupyter"
echo "c.ServerApp.allow_root = True" > "$VENV/etc/jupyter/jupyter_server_config.py"

if ! command -v myst >/dev/null 2>&1; then
  npm install -g mystmd --silent
fi

echo "export VIRTUAL_ENV=\"$VENV\"" >> "$CLAUDE_ENV_FILE"
echo "export PATH=\"$VENV/bin:\$PATH\"" >> "$CLAUDE_ENV_FILE"
