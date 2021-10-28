# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  run_Evolution.sh                                                                                         #
#  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                       #
#  6 Oct 2021                                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #    
#  Program to make and run the Evolution program and possibly reinstall the LHAPDF and/or Apfel++ library.  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  - The program also outputs the path of the currently used Apfel++ and LHAPDF Installation.               #
#    This output is aligned with he information output of the Evolution.cc file.                            #
#  - The variable $CURRENT_APFEL should point to the folder with the apfelxx folder inside.                 #
#  - The variable $CURRENT_LHAPDF should point to the folder with the LHAPDF-6.4.0 folder inside.           #
#  - The variable $$CURRENT_APFEL_TEST should point to the folder with the Evolution.cc and its             #
#    respectiveMakefile inside.                                                                             # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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


echo "Was the code in the apfelxx and/or LHAPDF folder changed?"
echo "Write 'both'/'b', 'apfel'/'a' or 'lhapdf'/'l' if so or anything else for no!"

changedCode=

while [[ $changedCode = "" ]]; do
    read changedCode
done

case $changedCode in
    both)
        installApfel
        installLHAPDF;;
    b)
        installApfel
        installLHAPDF;;
    apfel)
        installApfel;;
    a)
        installApfel;;
    lhapdf)
        installLHAPDF;;
    l)
        installLHAPDF;;
        
esac

echo ""
echo "Running the Evolution.cc program!"
echo ""

# Changing into directory
cd $CURRENT_APFEL_TEST

# temporary file to void unwanted console outputs 
voidOutput='voidOutput.txt' 

# Deleting previous Evolution program
make clean > $voidOutput

# Making new Evolution program
make Evolution > $voidOutput

# Running Evolution program
./Evolution

# Deleting previous Evolution program
make clean > $voidOutput

rm $voidOutput


echo "Used Apfel++ Installation : "$CURRENT_APFEL
echo "Used LHAPDF Installation  : "$CURRENT_LHAPDF

cd $CURRENT_APFEL && cd .. && cd plottingPython
python3 plottingEvolution.py &


echo ""
echo "Finished!!!"

cd $CURRENT_APFEL && cd .. 