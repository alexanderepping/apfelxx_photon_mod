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


# change, which programs should be installed again
changedCode="a"


case $changedCode in
    b)
        installApfel
        installLHAPDF;;
    a)
        installApfel;;
    l)
        installLHAPDF;;
esac

echo ""
echo "Running the StructureFunctions.cc program!"
echo ""

# Changing into directory
cd $CURRENT_APFEL_TEST

# temporary file to void unwanted console outputs 
voidOutput='voidOutput.txt' 

# Deleting previous StructureFunctions program
make clean > $voidOutput

# Making new StructureFunctions program
make StructureFunctions > $voidOutput

# Running StructureFunctions program
./StructureFunctions

# Deleting previous StructureFunctions program
make clean > $voidOutput

rm $voidOutput

cd $CURRENT_APFEL && cd .. && cd plottingPython
python3 plottingStructureFunctions.py &


echo ""
echo "Finished!!!"

cd $CURRENT_APFEL && cd .. 