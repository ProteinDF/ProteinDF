#!/bin/bash -x

BUILDDIR=build
SOURCEDIR=source

make gettext

# cp -R "${BUILDDIR}/locale/*.pot"  "${SOURCEDIR}/locale/pot"
# cp -R "${SOURCEDIR}/locale/pot/*" "${SOURCEDIR}/locale/ja/LC_MESSAGES/"
sphinx-intl update -p ${BUILDDIR}/locale -l en -l ja

