#!/bin/sh

git remote add -f tableimport2 git@github.com:glacier/tableimport2.git
git merge -s ours --no-commit tableimport2/master
git read-tree --prefix=tableimport -u tableimport2/master
