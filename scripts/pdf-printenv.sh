#!/bin/sh

if [ x${PDF_HOME} = x ]; then
    echo "unset PDF_HOME env. stop."
    exit 255
fi
echo "PDF_HOME=${PDF_HOME}"
