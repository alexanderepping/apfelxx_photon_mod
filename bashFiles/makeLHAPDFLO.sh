# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  makeLHAPDFLO.sh                                                                                    #
#  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                 #
#  22 Oct 2021                                                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Quick program to run the makeLHAPDFLO.f program and copy the file into the LHAPDF share directory. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cd ~/Uni/apfelxx_photon_mod/makeLHAPDF

gfortran -o makeLHAPDFLO makeLHAPDFLO.f grvphoton.f

./makeLHAPDFLO

rm ~/Uni/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSetLO/GRVCustomSetLO_0000.dat
rm ~/Uni/apfelxx_photon_mod/makeLHAPDF/makeLHAPDFLO
cp ~/Uni/apfelxx_photon_mod/makeLHAPDF/GRVCustomSetLO_0000.dat ~/Uni/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSetLO/
rm ~/Uni/apfelxx_photon_mod/makeLHAPDF/GRVCustomSetLO_0000.dat