# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  makeLHAPDFLO.sh                                                                                    #
#  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                 #
#  22 Oct 2021                                                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Quick program to run the makeLHAPDFLO.f program and copy the file into the LHAPDF share directory. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cd ~/Documents/apfelxx_photon_mod/makeLHAPDF

gfortran -o makeLHAPDF makeLHAPDFLO.f grvphoton.f

./makeLHAPDFLO

rm ~/Documents/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSetLO/GRVCustomSetLO_0000.dat
cp ~/Documents/apfelxx_photon_mod/makeLHAPDF/GRVCustomSetLO_0000.dat ~/Documents/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSetLO/