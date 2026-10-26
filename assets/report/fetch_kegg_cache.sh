#!/usr/bin/env bash
# Fill the rnaseq-splice report's KEGG cache with the KGML (.xml) and map image (.png) of every
# pathway in kegg_legacy_ids.tsv (MSigDB's 186 KEGG legacy gene sets), so reports can be rendered
# without internet access. Files already in the cache are skipped.
#
# Usage: fetch_kegg_cache.sh <cache_dir> [ids_tsv]
#   then render with: Rscript render_report.R ... --kegg_cache <cache_dir>
#
# KEGG's REST API is free for academic use; commercial use needs a KEGG licence (https://www.pathway.jp).
set -euo pipefail
cache=${1:?usage: fetch_kegg_cache.sh <cache_dir> [ids_tsv]}
ids=${2:-"$(dirname "$0")/kegg_legacy_ids.tsv"}
mkdir -p "$cache"
n=0; failed=0
for id in $(tail -n +2 "$ids" | cut -f2 | sort -u); do
    for kind in kgml:xml image:png; do
        what=${kind%%:*}; ext=${kind##*:}
        f="$cache/$id.$ext"
        [ -s "$f" ] && continue
        if curl -fsS --retry 2 -o "$f.part" "https://rest.kegg.jp/get/$id/$what"; then
            mv "$f.part" "$f"; n=$((n + 1))
        else
            rm -f "$f.part"; failed=$((failed + 1)); echo "failed: $id $what" >&2
        fi
        sleep 0.35   # KEGG asks for at most ~3 requests per second
    done
done
echo "Downloaded $n files into $cache ($failed failed)"
