#!/bin/bash
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=30:00:00
#PBS -l jobfs=400GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+gdata/xe2+gdata/fa63
#PBS -P xf3

set -xue

export PATH=/g/data/xf3/zz3507/app/RepeatMaskerV2/RepeatMasker_setup/RepeatMasker:/g/data/xf3/zz3507/app/RepeatMaskerV2/RepeatMasker_setup/RepeatModeler:$PATH

GENOME=/g/data/xf3/zz3507/Output/AfterHiCProcessing/2026Updatedgenome/H1/H1_2026_v1.top20.fa
LIB=/g/data/xf3/zz3507/Output/AfterHiCProcessing/H1V3H2V7ReferenceAUG17FFF/H1/RepeatModeler/123094750.gadi-pbs/AUG17H1V3FFF.500K.20scaffold-families.fa
OUTDIR=/g/data/xf3/zz3507/Output/AfterHiCProcessing/2026Updatedgenome/H1/RepeatMasker

mkdir -p "$OUTDIR"
cd "$OUTDIR"

RepeatMasker \
  -pa 48 \
  -e rmblast \
  -gff -xsmall \
  -lib "$LIB" \
  -dir "$OUTDIR" \
  "$GENOME"