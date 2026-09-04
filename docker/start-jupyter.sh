#!/bin/sh
set -eu

set -- jupyter lab \
  --ip=0.0.0.0 \
  --port="${JUPYTER_PORT:-8888}" \
  --no-browser \
  --allow-root \
  --ServerApp.root_dir="${JUPYTER_ROOT_DIR:-/workspace}"

if [ -n "${JUPYTER_TOKEN:-}" ]; then
  set -- "$@" --IdentityProvider.token="${JUPYTER_TOKEN}"
fi

exec "$@"
