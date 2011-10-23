#!/bin/sh

if [ `pwd` != "/Users/grace/Desktop/research" ]; then
    echo "You're not in the research directory .."
    exit 1
fi

if [ ! -e "/Volumes/research" ]; then
    echo "research from mbp is not mounted"
fi

rsync -av $1 --update . /Volumes/research
