#!/bin/bash
rm -rf bin lib
mkdir bin lib
cd bin
ln -s ../LinkDef_rdict.pcm
cd ..

cp Makefile.tmp Makefile

echo "now you can run \"make\" :)"
