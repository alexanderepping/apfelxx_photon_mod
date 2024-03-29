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
cd $CURRENT_APFEL_PROG 
cd ..
cd minimizationMinuit

echo ""
echo ""
echo "Compiling the Minimization program!"

# making the program
g++ -std=c++17 -I$CURRENT_APFEL/include -I$CURRENT_LHAPDF/include -I/usr/include/eigen3 -I/usr/include/boost -o minimizationStructureFunctions minimizationStructureFunctions.cc StructureFunctionsFcn.cc HessianMatrix.cc DeltaChi2.cc HelperFunctions.cc experimentalData.h -L/$CURRENT_APFEL/lib -lapfelxx -L/$CURRENT_LHAPDF/lib -lLHAPDF -L/usr/local/lib/ -lminuit-cpp 

echo ""
echo ""
echo "Running the Minimization program!"
echo ""

# Running Minimization program
./minimizationStructureFunctions


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