# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  makeLHAPDFHO.sh                                                                                    #
#  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                 #
#  25 Oct 2021                                                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Quick program to run the makeLHAPDFHO.f program and copy the file into the LHAPDF share directory. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cd ~/Uni/apfelxx_photon_mod/makeLHAPDF

gfortran -o makeLHAPDFHO makeLHAPDFHO.f grvphoton.f

./makeLHAPDFHO

rm ~/Uni/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSetHO/GRVCustomSetHO_0000.dat
rm ~/Uni/apfelxx_photon_mod/makeLHAPDF/makeLHAPDFHO
cp ~/Uni/apfelxx_photon_mod/makeLHAPDF/GRVCustomSetHO_0000.dat ~/Uni/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSetHO/
rm ~/Uni/apfelxx_photon_mod/makeLHAPDF/GRVCustomSetHO_0000.dat