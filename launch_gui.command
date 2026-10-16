#!/bin/bash
# Author: Khatereh Kashmari
# Date: 2026-08-10
# Version: 8.0.0

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FFBUILD_PYTHON="/usr/local/Caskroom/miniforge/base/envs/ffbuild/bin/python"
if [ ! -x "$FFBUILD_PYTHON" ]; then
  FFBUILD_PYTHON="/opt/homebrew/Caskroom/miniforge/base/envs/ffbuild/bin/python"
fi
if [ ! -x "$FFBUILD_PYTHON" ]; then
  echo "The ffbuild Python environment was not found."
  echo "Run parameterize_gui.py using a Python installation with Tkinter."
  read -r -p "Press Enter to close..."
  exit 1
fi
exec "$FFBUILD_PYTHON" "$SCRIPT_DIR/parameterize_gui.py"
