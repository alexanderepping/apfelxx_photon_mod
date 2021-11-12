
# To-Do:
[ ] write LaTeX file further
[ ] write README.md file for this folder

[ ] mainfile execution
[ ] `#pragma once` errors?
[ ] add timer and improve output
[ ] write Initial PDFs as function of the class StructureFunctionsFCN


[ ] config file: add parameters like (Output? y/n?)

[x] is it okay to leave the cutoff out? - Yes. The cutoff was only used because the output to the LHAPDF file couldn't be zero. If it doesn't produce any errors, we can just leave the negative numbers
[x] initial values? - Don't really matter
[x] add initial values in the config file


# Things to check, if there is some error:
- maybe introduce cutoff in StructureFunctionsFcn.cc in PDFsEvolved if there are some errors bc of negative numbers

# order of includes:
- pragma
- apfel, LHAPDF
- minuit
- myfiles
- std
	

# commands to compile etc:
see README.md