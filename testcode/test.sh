#!/bin/bash

trap "{echo quitting!; exit 1;}" SIGINT SIGTERM SIGKILL SIGHUP SIGSTOP
echo "pid is $$"

while :			# This is the same as "while true".
do
        sleep 60	# This script is not really doing anything.
done

