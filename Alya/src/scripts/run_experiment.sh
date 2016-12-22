#!/bin/bash

BASE="../.."
EXP_DIR="$BASE/experiments"
SRC_DIR="$BASE/src"

help_script()
{
    cat << EOF
Usage: $0 [options] /path_to_data_directory

Script for to get machine information before doing the experiment

OPTIONS:
   -h      Show this message
EOF
}

while getopts "d:" opt; do
    case $opt in
	d)
	    DATADIR="$OPTARG"
	    ;;
	h)
	    help_script
	    exit 4
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    help_script
	    exit 3
	    ;;
    esac
done

DATA_FOLD_DAY=`date +%Y_%m_%d`
DATA_FOLD_DAY="$EXP_DIR/$DATA_FOLD_DAY"
BKUP=`date +%H_%M_%S`
DATA_FOLD_HOST=`hostname`
DATA_FOLD_HOST="$DATA_FOLD_DAY/$DATA_FOLD_HOST"
DATA_FOLD_TIME="$DATA_FOLD_HOST/$BKUP"
INFO_FILE="$DATA_FOLD_TIME/Info.yaml"
DATA_FILE="$DATA_FOLD_TIME/Data.yaml"

mkdir -p $DATA_FOLD_DAY
mkdir -p $DATA_FOLD_HOST
mkdir -p $DATA_FOLD_TIME

ruby $SRC_DIR/boast/Split/run.rb --data=$DATA_FILE --info=$INFO_FILE --dimension=3
