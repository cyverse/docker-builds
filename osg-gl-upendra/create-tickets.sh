#!/usr/bin/env bash
  
set -e

mode="read"
while getopts ":rw" opt; do
    case "$opt" in
        r)
            mode="read"
            ;;
        w)
            mode="write"
            ;;
    esac
done
shift $((OPTIND -1))

path_file="$1"

echo "# application/vnd.de.path-list+csv; version=1"
for path in $(< "${path_file}"); do
    ticket=$(uuidgen | tr 'A-F' 'a-f' | tr -d '-' | fold -w 30 | head -n 1)
    iticket create "$mode" "$path" "$ticket"
    echo "${ticket},${path}"
done
