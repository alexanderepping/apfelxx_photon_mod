cd ~/Documents/apfelxx_photon_mod/makeLHAPDF
gfortran -o makeLHAPDF makeLHAPDF.f grvphoton.f
./makeLHAPDF
rm ~/Documents/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSet/GRVCustomSet_0000.dat
cp ~/Documents/apfelxx_photon_mod/makeLHAPDF/GRVCustomSet_0000.dat ~/Documents/apfelxx_photon_mod/myLHAPDFModified/share/LHAPDF/GRVCustomSet/