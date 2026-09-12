#!/bin/bash

#SCRIPT="/home/michaelgajda/thesis_work/isotope_fitting/my_flash_fitting_tppt_bin_optimization.py"
SCRIPT="C:/Users/burri/Documents/PET/Analysis_Programs/Analysis_Programs/Isotope_Fitting/tppt/my_flash_fitting_tppt_bin_optimization.py"
#PYTHON="/home/michaelgajda/mda/3d/.venv/bin/python"
PYTHON="C:/Users/burri/AppData/Local/Microsoft/WindowsApps/PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0/python"
#DATA_DIR="/home/michaelgajda/thesis_work/tppt_config"
DATA_DIR="C:/Users/burri/Documents/PET"
GROUP=4
DEBUG=false
if ${DEBUG};then

$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run6_TPPT_AcrylicNo2_Collimator3mm_30shots_15min_HWTrigOn_coinc.txt -group $GROUP
$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run7_TPPT_NoPhantom_Collimator2mm4mm8mm_30shots_15min_HWTrigOn_coinc.txt -group $GROUP
$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run8_TPPT_AcrylicNo3_Collimator2mm4mm8mm_30shots_15min_HWTrigOn_coinc.txt -group $GROUP
$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run9_TPPT_Acrylic3mm_NoCollimator_30shots_15min_HWTrigOn_coinc.txt -group $GROUP
$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run10_TPPT_Acrylic2mm4mm8mm_NoCollimator_30shots_15min_HWTrigOn_coinc.txt -group $GROUP
$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run11_TPPT_AcrylicNo4_Delta6Collimator_30shots_15min_HWTrigOn_coinc.txt -group $GROUP
$PYTHON $SCRIPT -f $DATA_DIR/Run4_TPPT_AcrylicEmptyPipe.txt -group $GROUP

$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run2_FullAcrylicNo2_Collimator19.75cm_89.6MeV_HWTrigOn_15min_coinc.txt -group $GROUP

$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run3_FullAcrylicNo3_Collimator19.75cm_103.8MeV_HWTrigOn_15min_coinc.txt -group $GROUP
$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run5_TPPT_Acrylic1_NoCollimator_1shot_15min_HWTrigOn_coinc.txt -group $GROUP
else

#$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run1_FullAcrylicNo1_Collimator19.75cm_72.5MeV_HWTrigOn_15min_coinc.txt -group $GROUP

#$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run2_FullAcrylicNo2_Collimator19.75cm_89.6MeV_HWTrigOn_15min_coinc.txt -group $GROUP

#$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run3_FullAcrylicNo3_Collimator19.75cm_103.8MeV_HWTrigOn_15min_coinc.txt -group $GROUP
#$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run5_TPPT_Acrylic1_NoCollimator_1shot_15min_HWTrigOn_coinc.txt -group $GROUP

$PYTHON $SCRIPT -f $DATA_DIR/FilteredData_Run6_Full1inGraphite_Collimator19.75cm_103.8MeV_HWTrigOn_15min_coinc.txt -group $GROUP

fi

