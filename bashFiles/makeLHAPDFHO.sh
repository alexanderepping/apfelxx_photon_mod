# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  makeLHAPDFHO.sh                                                                                    #
#  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                 #
#  25 Oct 2021                                                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Quick program to run the makeLHAPDFHO.f program and copy the file into the LHAPDF share directory. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cd ~/Documents/apfelxx_photon_mod/makeLHAPDF

gfortran -o makeLHAPDFHO makeLHAPDFHO.f grvphoton.f

./makeLHAPDFHO

rm ~/Documents/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSetHO/GRVCustomSetHO_0000.dat
rm ~/Documents/apfelxx_photon_mod/makeLHAPDF/makeLHAPDFHO
cp ~/Documents/apfelxx_photon_mod/makeLHAPDF/GRVCustomSetHO_0000.dat ~/Documents/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSetHO/
rm ~/Documents/apfelxx_photon_mod/makeLHAPDF/GRVCustomSetHO_0000.dat