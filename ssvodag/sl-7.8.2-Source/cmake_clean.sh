#!/bin/sh

make clean
find . -name "install_manifest.txt" -exec rm -fr {} \;
find . -name "CMakeCache.txt" -exec rm -fr {} \;
find . -name "*.cmake" -exec rm -fr {} \;
find . -name "CMakeFiles" -exec rm -fr {} \;
find . -name "_CPack_Packages" -exec rm -fr {} \;
find . -name "Makefile" -exec rm -fr {} \;
find . -name "Testing" -exec rm -fr {} \;
