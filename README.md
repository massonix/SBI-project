# SBI-Python project  
## Constructing macromolecular complexes. 

*Ramon Massoni and Winona Oliveros*

### Description
A python package for macrocomplex construction given protein pairwise interactions. 

**Problem**

Given a set of interacting pairs (prot-prot), reconstruct the complete macrocomplex and return a PDB file (or files) with the possible protein macrocomplexes build. 

## General Information
### Input Files

This program needs an input of PDB files holding the protein pairwise interactions needed to reconstruct the desired macrocomplex. The program can handle those scenarios: 

* The same sequence appearing in different PDB files has not to be identical, we can handle 95% of identity. 
* The same sequence appearing in different files can have different names. 

##Background and algorithm

Many proteins consists of several polypeptide chains that assembly into multi-subunit complexes. Such assembly is known
as the quaternary structure of a protein, and it plays a major role in defining a protein's function, localization and activity. Examples include proteins as important as the nucleosome, the spliceosome or the proteasome, among others. Nevertheless, understanding how individual subunits assemble into larger complexes is no easy task. In this project, we developed a stand-alone application to create multi-subunit complexes out of its individual pairwise interactions. 

We approached the former problem by strucurally superimposing structures that had at least two identical subunits (polypeptide chains that have a pairwise sequence identity >= 95%). We considered as feasible complexes those superimpositions that had no clashes, which means that the backbone of one structure did not protrude the backbone of the other. We then called the function again giving as input the output of the first call (i.e. the resulting complex), establishing in this manner a recursive algorithm. Nevertheless, as several complexes could be found from the same input, the computational-cost grew exponentially. To overcome such challenge, we applied two restraints.

First, as the same complex could be achieved following different paths, we ensured that each new complex was unique before calling the function again. Again, we did this by superimposing each complex with all the former ones. If the complex was repeated, at least one superimposition resulted in an extremely low RMSD (< 2). In such cases, the complex was discarded.

Second, we were very restrictive with the stoichiometry of the complex. The stoichiometry of a protein complex refers to the composition of its subunits, and it is specified by a formula. Such formula consists of uppercase letters with integer coefficients. The letters refer to unique subunits, and the coefficients to the amount of each unique subunit. For instance, hemoglobin has two alpha and two beta subunits, and thus a stoichiometry of A2B2. Larger coeffiecients must be used first, i.e A2B3 is not valid, A3B2 is. We ask the user to provide the stoichiometry he or she desires to obtain in the final complex, and hence all the complexes that do not match such stoichiometry are ruled out, achieving a more efficient program. Only those complexes which meet the stoichiometry are saved in individual PDB files. 

Finally, we also allow the user to optimize the resulting complexes, which means that the side chain and loops may be restructured to minimize the energy and end up with a more reliable complex.


##Tutorial

The following section explains the user how to use our application to make the most out of it.

###Command-line arguments

* -i --input: directory containing the PDB files that will serve as input to the program.                        
* -o --output: directory where the macrocomplexes will be saved as PDB files
* -v --verbose: whether or not the user wants to keep track of what the programme is doing.                       
* -st --stoich: the stoichiometry of the multi-subunit complexes.                     
* -opt --optimization: if the user wants to optimize the resulting complexes.        
* -cn --contact_num: number amino acids that are allowed to protrude after each superimposition (default=5).
* -f --files_num: maximum number of files (complexes) the user wants to obtain.                   

To get a better understanding of how to run the programme properly, we show 3 different scenarios that represent the heterogeneity of inputs that may be provided.

####Example 1

The first example corresponds to the phosphate dehydratase (pdb id: 2F1D). As we can see in the splited file we provide, (2f1d_split.pdb), this is an Homo 8-mer (stoichiometry: A8). To achieve this complex we can run:

```python
python3 multicomplex.py -i example1 -o complexes/ -st "A8" -f10 -v
```
Where example1 is the directory containing all input files, complexes is the directory where the output files will be saved, "A8" is the stoichiometry (i.e. 8 equal chains), -f10 means that the programme will stop after 10 complexes, and -v means that the standard error will be printed. In this case, we can see how multiple complexes are output. This is because at the recurssion level the function stopped (k=8), the programme could still superimpose more structures without finding clashes. Interestingly, we see how the complex9.pdb is the same as the 2f1d_split (inside raw_pdbs), except for one chain.

As more chains can be added, we can rerun the programme with a higher stoichiometry: A24. This time though, it will take longer to execute, so we will limit the max number of files to 1. Before doing so, we must move all the previous complexes to a different file or change the outdir:

```python
python3 multicomplex.py -i example1 -o complexes/ -st "A24" -f1 -v
```
As in this example we only have one unique chain, if we try to get a hetero k-mer, we will get an error:

```python
python3 multicomplex.py -i example1 -o complexes/ -st "A2B2" -f1 -v
```

####Example 2

The second example was kindly provided by Prof. Baldo Oliva. In this case, we have two unique chains, so we can build hetero k-mers:

```python
python3 multicomplex.py -i example2 -o complexes -st "A6B6" -f1 -v 
```
As we can see, it completes the recursion pretty fast. However, if we try to add two more chains the resulting complexes will clash, so it will take too long and do not give any result:

```python
python3 multicomplex.py -i example2 -o complexes -st "A7B7" -f1 -v 
```
Thus, if we are not interested in any particular stoichiometry, it is a good practice to start with few recursions, check the output, and increase it or not accordingly. Finally, in the case were we have more than one unique chain, we can still try to find Homo k-mers

```python
python3 multicomplex.py -i example2 -o complexes -st "A6" -f1 -v 
```

###Example 3

The third example corresponds to the 20S proteasome (pdb id: 1G65). This represents an extreme case, where every chain is different. The proteasome is a Hetero 28-mer (stoichiometry: A2B2C2D2E2F2G2H2I2J2K2L2M2N2). For our purposes, we splited in two and obtained a stoichiometry of: A1B1C1D1E1F1G1H1I1J1K1L1M1N1 (1g65_split, inside raw_pdbs). We will input increasing stoichiometries and assess the scalability of the programme accordingly:

```python
python3 multicomplex.py -i example3 -o complexes -st "A1B1C1D1" -f4 -v 
python3 multicomplex.py -i example3 -o complexes -st "A1B1C1D1E1F1G1" -f4 -v 
python3 multicomplex.py -i example3 -o complexes -st "A1B1C1D1E1F1G1H1I1J1K1L1M1N1" -f4 -v 
```
The last example gives exactly half 20S proteasome. However, it took extremely long to execute (~15min). Thus, if one wants to obtain a complex with a great deal of unique chains, it is a wise idea to let the programme run and go have a coffee ;).

##Limitations

* Unable to model protein-RNA or protein-DNA interactions.
* Scalability: after certain recursion depth the program does not scalate too well. To improve the program, we could have used a dynammic programming approach, where out of several complexes only the most optimal is kept, and the next solutions are built onto that one.
* Stoichiometry: as powerful as it may be, the application relies too heavily on the stoichiometry input. On the one hand, this gives a greater control to the user over the output he or she expects but, on the other hand, it does not find complexes blindly, regardless of the stoichiometry.







