#!/bin/bash

for path in */; do
    [ -d "${path}" ] || continue # if not a directory, skip
    dirname="$(basename "${path}")"

    cd $dirname

    make
    sudo make install
    make clean

    cd ..
done