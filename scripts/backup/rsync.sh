#!/bin/bash

rsync -av --progress --update $1 /media/disk/work_new/week_* .
