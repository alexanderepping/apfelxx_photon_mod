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

echo ""
echo ""
echo "Compiling the Minimization program!"

# making the program
time g++ -std=c++17 -I/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/include -I/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/include -I/usr/include/eigen3 -I/usr/include/boost -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc ErrorPDFs.cc OutputFunctions.cc ResultsFunctions.h experimentalData.h -L/home/alexander/Uni/apfelxx_photon_mod/myApfelModified/lib -lapfelxx -L/home/alexander/Uni/apfelxx_photon_mod/myLHAPDFModified/lib -lLHAPDF -L/home/alexander/Uni/minuit-cpp/build/lib/ -lminuit-cpp 

echo ""
echo ""
echo "Running the Minimization program!"
echo ""

# Running Minimization program
time ./minimizationStructureFunctions


# deleting the program
rm ./minimizationStructureFunctions

echo "Used Apfel++ Installation : "$CURRENT_APFEL
echo "Used LHAPDF Installation  : "$CURRENT_LHAPDF

echo "Plot Initial PDFs, including Errors? (y or n)"
echo "Depending on the In-/OutputFiles of the different programs, the plot might not be the calculated one."

answer=""

while [[ $answer = "" ]]; do
    read answer
done

echo ""

if [ $answer = "y" ]; then
    cd $CURRENT_APFEL && cd .. && cd plottingPython
    python3 plottingInitialPDFs.py &
fi


echo "Calculate the Structure Functions, including Errors? (y or n)"
echo "Depending on the In-/OutputFiles of the different programs, the result might not be the calculated one."

answer=""

while [[ $answer = "" ]]; do
    read answer
done

echo ""

if [ $answer = "y" ]; then
    cd $CURRENT_APFEL && cd .. && cd bashFiles
    bash run_StructureFunctions.sh
fi


echo ""
echo "Finished!!!"

cd $CURRENT_APFEL && cd .. 