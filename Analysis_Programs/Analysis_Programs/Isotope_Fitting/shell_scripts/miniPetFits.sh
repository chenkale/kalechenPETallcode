#!/bin/bash

SCRIPT="/home/michaelgajda/thesis_work/isotope_fitting/my_FLASH_Fitting.py"
PYTHON="/home/michaelgajda/mda/3d/.venv/bin/python"
DATA_DIR="/home/michaelgajda/thesis_work/mini_isotope_config"
GROUP=4
DEBUG=true

if ${DEBUG}; then
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-Cu-10min-210mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-H2O1-15min-210mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-HDPE1-20min-210mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-Ni-10min-210mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-NoPhantom-5min-210mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-PMMA1-3min-210mm-4cmCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-PMMA1-3min-210mm-8cmCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-PMMA1-3min-210mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-PMMA2-20min-210mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-PMMA3-15min-150mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-PMMA3-15min-150mm-NoCu-2.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-PMMA4-15min-110mm-NoCu.txt -group $GROUP
  $PYTHON $SCRIPT -f $DATA_DIR/MDA-PMMA5-20min-210mm-NoCu.txt -group $GROUP
else
  # You can add your own else clause if needed, or leave it empty
  echo "No non-debug files specified."
fi

