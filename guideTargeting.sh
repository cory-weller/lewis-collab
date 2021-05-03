#!/usr/bin/env bash

module load python/3.6

python ./guideTargeting.py ${1} ${2} ${3} > data/guideTargetCounts.tsv