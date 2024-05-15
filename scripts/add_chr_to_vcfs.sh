#!/usr/bin/env bash

# This bash script adds a 'chr' to the chromosome names in a vcf.bgz file, then bgzips the output.

zcat $1 | \
 awk '{
        if($0 !~ /^#/)
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0
      }' | \
 bgzip -c > $2