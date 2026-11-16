#!/usr/bin/env bash
# Put the built Sphinx HTML into a gh-pages checkout, under python-docs/ only.
#
#   tools/publish_docs.sh <built html dir> <gh-pages checkout>
#
# The gh-pages branch holds the DISCUS homepage. This script replaces the
# python-docs/ folder and makes sure an empty .nojekyll exists at the root
# (GitHub Pages' Jekyll would otherwise drop Sphinx's _static/, _images/, ...;
# the homepage is plain HTML and is not affected). It then checks that
# nothing else in the checkout changed, and stages the change. Committing and
# pushing are left to the caller (CI).
#
# Exit codes: 0 = staged (or nothing to do), 1 = usage error,
#             2 = something outside python-docs/ and .nojekyll would change.
set -euo pipefail

DEST=python-docs

if [ $# -ne 2 ] || [ ! -f "$1/index.html" ] || ! git -C "$2" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "usage: $0 <built html dir with index.html> <gh-pages git checkout>" >&2
    exit 1
fi
html=$(cd "$1" && pwd)
site=$(cd "$2" && pwd)

cd "$site"
rm -rf -- "$DEST"
mkdir -p -- "$DEST"
cp -a -- "$html"/. "$DEST"/
rm -f -- "$DEST/.buildinfo"          # Sphinx build cache, not needed online
[ -e .nojekyll ] || : > .nojekyll

git add -A -- "$DEST" .nojekyll

# Guard: only python-docs/ and .nojekyll may differ from the last commit.
outside=$(git status --porcelain | awk '{print $NF}' \
          | grep -v -e "^$DEST/" -e "^$DEST\$" -e '^\.nojekyll$' || true)
if [ -n "$outside" ]; then
    echo "refusing to publish: these paths outside $DEST/ would change:" >&2
    echo "$outside" >&2
    exit 2
fi

if git diff --cached --quiet; then
    echo "docs unchanged; nothing to publish"
else
    git diff --cached --stat | tail -1
fi
