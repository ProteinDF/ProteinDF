#!/bin/bash


build_html()
{
    LANG=$1
    if [ x${LANG} = x ]; then
        LANG="en"
    fi
    cd source
    sphinx-intl build -l "${LANG}"
    sphinx-build -b html -a -E -D language="${LANG}" . ../build/html/${LANG}
}


echo ">>>> build en"
build_html en

echo ">>>> build ja"
build_html ja
