#!/bin/bash
# memguard.sh LIMIT_MB SECONDS cmd... : runs cmd, kills its process tree if RSS > LIMIT_MB or time > SECONDS.
limit=$1; secs=$2; shift 2
"$@" &
pid=$!
start=$(date +%s)
tree() { local p=$1; echo $p; for c in $(pgrep -P $p); do tree $c; done; }
while kill -0 $pid 2>/dev/null; do
  rss=0; for p in $(tree $pid); do r=$(ps -o rss= -p $p 2>/dev/null); rss=$((rss + ${r:-0})); done
  if [ $((rss / 1024)) -gt $limit ] || [ $(( $(date +%s) - start )) -gt $secs ]; then
    echo "MEMGUARD: killing (rss=$((rss/1024))MB, t=$(( $(date +%s) - start ))s)" >&2
    for p in $(tree $pid); do kill -9 $p 2>/dev/null; done
    wait $pid 2>/dev/null; exit 137
  fi
  sleep 0.3
done
wait $pid
