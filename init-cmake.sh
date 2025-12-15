#!/bin/bash
rm -rf build bin
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
