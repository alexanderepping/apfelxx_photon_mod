
# To-Do:
[ ] more data points

[ ] makefile or sth like that


[x] is it okay to leave the cutoff out? - Yes. The cutoff was only used because the output to the LHAPDF file couldn't be zero. If it doesn't produce any errors, we can just leave the negative numbers
[x] initial values? - Don't really matter
[x] add initial values in the config file
[x] write Initial PDFs as function of the class StructureFunctionsFCN
[x] mainfile execution
[x] change README.md file for this folder 
[x] improve output
[x] momentum sum rule by itself
[x] define the parameters inside the initial pdfs before so that I can always use the same functions
[x] option for IntitalPDFs to print or not print
[x] switch case in main for parameters
[x] file output, such that I can plot it with python
[x] python file to plot:
[x] write LaTeX file further
[x] `#pragma once` errors? - Not important. Deleting them leads to errors.


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