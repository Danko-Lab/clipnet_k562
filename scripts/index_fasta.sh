#/usr/bin/env bash

# Adds record number to each header
awk '{if (/^>/) print ">"(++i)"_" substr($0,2); else print $0;}' $1

# Replaces header with record number
#awk '/^>/{print ">" ++i; next}{print}' $1