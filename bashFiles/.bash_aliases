#std, mod
apfelInstallation="mod"
lhapdfInstallation="mod"



###################
# apfelxx
###################
case $apfelInstallation in
"std")
export LD_LIBRARY_PATH=/local0/a_eppi01/myapfel/lib:$LD_LIBRARY_PATH
export PATH=/local0/a_eppi01/myapfel/bin:$PATH
export CURRENT_APFEL=/local0/a_eppi01/myapfel
export CURRENT_APFEL_TEST=/local0/a_eppi01/test_calculations ;;

"mod")
export LD_LIBRARY_PATH=/local0/a_eppi01/myApfelModified/lib:$LD_LIBRARY_PATH
export PATH=/local0/a_eppi01/myApfelModified/bin:$PATH
export CURRENT_APFEL=/local0/a_eppi01/myApfelModified
export CURRENT_APFEL_TEST=/local0/a_eppi01/testModifiedApfel;;
esac



###################
# lhapdf
###################
case $lhapdfInstallation in
"std")
export LD_LIBRARY_PATH=/local0/a_eppi01/myLHAPDF/lib:$LD_LIBRARY_PATH
export PATH=/local0/a_eppi01/myLHAPDF/bin:$PATH
export CURRENT_LHAPDF=/local0/a_eppi01/myLHAPDF ;;

"mod")
export LD_LIBRARY_PATH=/local0/a_eppi01/myLHAPDFModified/lib:$LD_LIBRARY_PATH
export PATH=/local0/a_eppi01/myLHAPDFModified/bin:$PATH
export CURRENT_LHAPDF=/local0/a_eppi01/myLHAPDFModified;;
esac
# maybe needed:
# export PYTHONPATH=$PYTHONPATH:/foo/lhapdf/lib/python3.9/site-packages
