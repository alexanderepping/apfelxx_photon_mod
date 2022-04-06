
# To-Do:
[x] look through the mathematica docs wo see which are needed
[ ] make the descriptions in pointlike... correct and also include the sources
[x] git commit
[ ] implement, that I don't have to manually change the pointlike perturbation order


[ ] what about my apfel NLO results? compare them to LO results (apfel and GRV) & compare the GRV LO and NLO results
[ ] what could be changed so that I get more similar values to the ones given by SAL? 
[ ] same results as Vadim
[ ] to check how good the values for the parameters are we can look at the chi2 for the different experiments seperately

[ ] makefile or sth like that

[ ] I won't need to do the dijet stuff, I can just treat it as a black box an d vadim will send me the stuff I have to do


## done To-Dos:
[x] implement momentum sum rule FG (pointlike part is left)
[x] run the minimization w/ the SAL PDFs 
        - mainly test if the new momentum sum rule is working
[x] add bounds to SAL initial pdfs


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

# for commits:
- there was also an error w/ the previous momentumsumruleSAL, because params[0] should've been counted twice

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
