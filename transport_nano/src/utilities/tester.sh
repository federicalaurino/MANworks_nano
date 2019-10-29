#!/bin/bash

inputfile=input_Simulation_A.param
run="./M3D1D $inputfile"

echo "========================================"
echo "RIGID STRONG 0.9"
echo "========================================"
sed -r -i "s/^RIGID_STRONG.*/RIGID_STRONG = 1;/" $inputfile
sed -r -i "s/^rho_l.*/rho_l = 0.9;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/rigid_strong_09
printfile=rigid_strong_09.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
sed -r -i "s/^RIGID_STRONG.*/RIGID_STRONG = 0;/" $inputfile

echo "========================================"
echo "RIGID MILD 0.3"
echo "========================================"
sed -r -i "s/^RIGID_MILD.*/RIGID_MILD = 1;/" $inputfile
sed -r -i "s/^rho_l.*/rho_l = 0.3;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/rigid_mild_03
printfile=rigid_mild_03.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "RIGID MILD 0.7"
echo "========================================"
sed -r -i "s/^rho_l.*/rho_l = 0.7;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/rigid_mild_07
printfile=rigid_mild_07.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "RIGID MILD 0.9"
echo "========================================"
sed -r -i "s/^rho_l.*/rho_l = 0.9;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/rigid_mild_09
printfile=rigid_mild_09.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

sed -r -i "s/^RIGID_MILD.*/RIGID_MILD = 0;/" $inputfile


echo "========================================"
echo "SOFT MILD 0.3"
echo "========================================"
sed -r -i "s/^SOFT_MILD.*/SOFT_MILD = 1;/" $inputfile
sed -r -i "s/^rho_l.*/rho_l = 0.3;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/soft_mild_03
printfile=soft_mild_03.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "SOFT MILD 0.7"
echo "========================================"
sed -r -i "s/^rho_l.*/rho_l = 0.7;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/soft_mild_07
printfile=soft_mild_07.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "SOFT MILD 0.9"
echo "========================================"
sed -r -i "s/^rho_l.*/rho_l = 0.9;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/soft_mild_09
printfile=soft_mild_09.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
sed -r -i "s/^SOFT_MILD.*/SOFT_MILD = 0;/" $inputfile

echo "========================================"
echo "SOFT STRONG 0.3"
echo "========================================"
sed -r -i "s/^SOFT_STRONG.*/SOFT_STRONG = 1;/" $inputfile
sed -r -i "s/^rho_l.*/rho_l = 0.3;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/soft_strong_03
printfile=soft_strong_03.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "SOFT MILD 0.7"
echo "========================================"
sed -r -i "s/^rho_l.*/rho_l = 0.7;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/soft_strong_07
printfile=soft_strong_07.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk

echo "========================================"
echo "SOFT MILD 0.9"
echo "========================================"
sed -r -i "s/^rho_l.*/rho_l = 0.9;/" $inputfile
resu_folder=./pressure_5mmHg/cinj_0_1/soft_strong_09
printfile=soft_strong_09.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mkdir vtk
sed -r -i "s/^SOFT_STRONG.*/SOFT_STRONG = 0;/" $inputfile


