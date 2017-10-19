#!/bin/bash

BASE="."
EXP_DIR="$BASE/experiments"
SRC_DIR="$BASE"

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
KERNEL_FILE="$DATA_FOLD_TIME/kernel.f90"

mkdir -p $DATA_FOLD_DAY
mkdir -p $DATA_FOLD_HOST
mkdir -p $DATA_FOLD_TIME

#echo "EXECUTING ruby $SRC_DIR/run.rb --data=$DATA_FILE --kernel=$KERNEL_FILE"
#ruby $SRC_DIR/run.rb --data=$DATA_FILE --kernel=$KERNEL_FILE


echo "EXECUTING ruby $SRC_DIR/run.rb --data=$DATA_FILE --kernel=$KERNEL_FILE"
ruby $SRC_DIR/run_compOPT.rb --data=$DATA_FILE --kernel=$KERNEL_FILE
echo "EXECUTING ruby parse_data_yaml.rb $DATA_FILE test.cvs"
ruby parse_data_yaml.rb $DATA_FILE test.cvs
echo "EXECUTING sh parse_cvs.sh"
sh parse_cvs.sh
