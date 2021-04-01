# Density Matrix Embedding Theory (DMET)
DMET and one-shot DMET codes: CISD-in-HF (high-level and low-level calculations)

Copyright &copy; 2016 Borna Zandkarimi  

## **General comments regarding the code**:

1. All calculations are done in atomic units.  
2. Contracted Gaussian functions are used as basis sets.  
3. The basis set should be in a file named "Basis.txt".  
4. You just need to change the name of prepared basis set files to “Basis.txt”, or copy the text of these prepared basis set files to the “Basis.txt”.  
5. The geometry will be saved in a file named "Input.xyz" in the XYZ format.  
6. The final results will be saved in a file named "Output.txt".  
7. The initial values of input data can be changed in the “allocate-array.f90”.  
8. LAPACK library needs to be installed. The path to it can be changed in the "make" file.  
