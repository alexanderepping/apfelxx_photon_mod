# Changing into directory
cd $CURRENT_APFEL_TEST

# temporary file to void unwanted console outputs 
voidOutput='voidOutput.txt'

# Deleting previous programs
make clean > $voidOutput

# Making new StructureFunctionsErrorPDFs program
make StructureFunctionsErrorPDFs > $voidOutput


# Running StructureFunctionsErrorPDFs program
./StructureFunctionsErrorPDFs

make clean > $voidOutput
rm $voidOutput

cd $CURRENT_APFEL && cd .. && cd plottingPython
python3 plottingStructureFunctionsErrorPDFs.py &

cd $CURRENT_APFEL && cd .. 