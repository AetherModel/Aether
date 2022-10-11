#!/bin/sh

# Usage:
#     ./run.sh {example}
# Example:
#     ./run.sh sphere_2d

mpic++ -I /usr/local/include \
       -I ../../../share/include \
       -I ../../../include \
       -o "$1" "$1".cpp
