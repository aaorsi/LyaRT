#!/bin/bash
HOST=$1
shift
for ARG in "$@"
do
        nmap -p $ARG $HOST
done
