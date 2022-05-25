# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  run_EvolutionFlavors.sh                                                                                          #
#  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                               #
#  22 Oct 2021                                                                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #    
#  Program to make and run the minimizationMinuit program and possibly reinstall the LHAPDF and/or Apfel++ library. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  - The program also outputs the path of the currently used Apfel++ and LHAPDF Installation.                       #
#    This output is aligned with he information output of the EvolutionFlavors.cc file.                             #
#  - The variable $CURRENT_APFEL should point to the folder with the apfelxx folder inside.                         #
#  - The variable $CURRENT_LHAPDF should point to the folder with the LHAPDF-6.4.0 folder inside.                   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

installApfel () {    
    #running the commands to install the package again
    cd $CURRENT_APFEL/apfelxx/build
    cmake -DCMAKE_INSTALL_PREFIX=$CURRENT_APFEL/ ..
    make && make install
    
    echo ""
    echo "The code for Apfel++ was installed again!"
}

installLHAPDF () {    
    #running the commands to install the package again
    cd $CURRENT_LHAPDF/LHAPDF-6.4.0
    ./configure --prefix=$CURRENT_LHAPDF/ 
    make && make install

    echo ""
    echo "The code for LHAPDF was installed again!"
}

installApfel

# Changing into directory
cd $CURRENT_APFEL_TEST 
cd ..
cd minimizationMinuit

# Deleting previous program
rm ./testMinimization

echo ""
echo ""
echo "Compiling the Minimization program!"

# making the program
time g++ -std=c++17 -I/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/include -I/usr/include/eigen3 -o testMinimization minimizationStructureFunctions.cc StructureFunctionsFcn.cc experimentalData.h ErrorPDFs.cc -L/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF -L/home/alexander/Uni/minuit-cpp/build/lib/ -lminuit-cpp 

echo ""
echo ""
echo "Running the Minimization program!"
echo ""

# Running Minimization program
time ./testMinimization

echo "Used Apfel++ Installation : "$CURRENT_APFEL
echo "Used LHAPDF Installation  : "$CURRENT_LHAPDF

cd $CURRENT_APFEL && cd .. && cd plottingPython
python3 plottingInitialPDFsWithError.py &


echo ""
echo "Finished!!!"

cd $CURRENT_APFEL && cd .. 