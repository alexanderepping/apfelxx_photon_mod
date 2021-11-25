
# To-Do:
[ ] makefile or sth like that

[ ] what about my apfel NLO results? compare them to LO results (apfel and GRV) & compare the GRV LO and NLO results

[ ] implement momentum sum rule FG (pointlike part is left)
[ ] add bounds to SAL initial pdfs


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
[x] more data points
[x] take a look at vadims notes
[x] take a look at SAL paper
[x] take a look at nCTEQ15 section 2 for error analysis
[x] more experimental data
[x] correct the energies (used Q2 instead of Q)
[x] correct the momentum sum rule 
[x] add all the stuff with the enumerator etc
[x] implement the SAL inputPDFs
[x] add initialPDFs from SAL to python


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