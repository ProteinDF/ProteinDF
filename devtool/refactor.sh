#!/bin/bash

set -x

from_str="${1}"
to_str="${2}"


usage()
{
    echo "${0} <from_str> <to_str>"
}

if [ "x${from_str}" = x -o "x${to_str}" = x ]; then 
    usage
    exit
fi

echo "from_str=${from_str}"
echo "to_str=${to_str}"

find . \( -name '*.h' -or -name '*.cpp' -or -name '*.cc' -or -name 'CMakeLists.txt' \) \
  -print \
  -exec sed -i -e "s/${from_str}/${to_str}/g" {} \;
