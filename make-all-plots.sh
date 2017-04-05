#!/usr/bin/env bash

if [[ $- == *i* ]]; then
    echo "don't source" >&2
    return 1
else
    set -eu
fi

EXT='.pdf'

function _usage() {
    echo "usage: ${0##*/} [-h] [opts]"
}

function _help() {
    _usage
    cat <<EOF
Abandon all hope
 -h: print help
 -s: draw slices
 -e <ext>: plot extension (default ${EXT})
EOF
}

DRAW_SLICES=''
while getopts ":he:s" opt $@; do
    case $opt in
        h) _help; exit 1;;
        e) EXT=${OPTARG} ;;
        s) DRAW_SLICES=1 ;;
        # handle errors
        \?) _usage; echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :) _usage; echo "Missing argument for -$OPTARG" >&2; exit 1;;
        *) _usage; echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done

echo "drawing exclusion contour"
draw-limits.py sigma_theory.dat sigma_exp.dat sigma_sigma_up1.dat sigma_sigma_do1.dat -e $EXT

if [[ ! $DRAW_SLICES ]]; then
    exit 0
fi

POINTS=$(cat sigma_theory.dat  | awk '{print $1}' | sort -ug)

for POINT in ${POINTS}; do
    echo "making slice $POINT"
    draw-slice.py sigma_theory.dat sigma_exp.dat -z $POINT -e $EXT
done



