#!/bin/bash

for extension in "*.cpp" "*.c" "*.h" "*.hpp"; do
    find source -name $extension -exec \
        clang-format {} -i -style="{IndentWidth: 4, ColumnLimit: 79}" \;
done
