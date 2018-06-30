#!/bin/bash

find . \( -name '*.h' -or -name  '*.cc' -or -name '*.cpp' \) -type f -exec clang-format -i {} \;


