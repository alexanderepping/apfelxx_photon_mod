#std, mod
apfelInstallation="mod"
lhapdfInstallation="mod"



###################
# apfelxx
###################
case $apfelInstallation in
"std")
export LD_LIBRARY_PATH=/local0/a_eppi01/apfelxx_photon_mod/myapfel/lib:$LD_LIBRARY_PATH
export PATH=/local0/a_eppi01/apfelxx_photon_mod/myapfel/bin:$PATH
export CURRENT_APFEL=/local0/a_eppi01/apfelxx_photon_mod/myapfel
export CURRENT_APFEL_TEST=/local0/a_eppi01/apfelxx_photon_mod/test_calculations ;;

"mod")
export LD_LIBRARY_PATH=/local0/a_eppi01/apfelxx_photon_mod/myApfelModified/lib:$LD_LIBRARY_PATH
export PATH=/local0/a_eppi01/apfelxx_photon_mod/myApfelModified/bin:$PATH
export CURRENT_APFEL=/local0/a_eppi01/apfelxx_photon_mod/myApfelModified
export CURRENT_APFEL_TEST=/local0/a_eppi01/apfelxx_photon_mod/testModifiedApfel;;
esac



###################
# lhapdf
###################
case $lhapdfInstallation in
"std")
export LD_LIBRARY_PATH=/local0/a_eppi01/apfelxx_photon_mod/myLHAPDF/lib:$LD_LIBRARY_PATH
export PATH=/local0/a_eppi01/apfelxx_photon_mod/myLHAPDF/bin:$PATH
export CURRENT_LHAPDF=/local0/a_eppi01/apfelxx_photon_mod/myLHAPDF ;;

"mod")
export LD_LIBRARY_PATH=/local0/a_eppi01/apfelxx_photon_mod/myLHAPDFModified/lib:$LD_LIBRARY_PATH
export PATH=/local0/a_eppi01/apfelxx_photon_mod/myLHAPDFModified/bin:$PATH
export CURRENT_LHAPDF=/local0/a_eppi01/apfelxx_photon_mod/myLHAPDFModified;;
esac
# maybe needed:
# export PYTHONPATH=$PYTHONPATH:/foo/lhapdf/lib/python3.9/site-packages
