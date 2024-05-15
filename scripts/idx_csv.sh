#!/bin/bash

awk '{$0=(++i)"_"$0}1' $1