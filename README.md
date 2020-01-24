# SPECS

SPECS is a side-chain-orientation-included protein model-native similarity metric for improved evaluation of protein structural models. SPECS stands for <b>Superposition-based Protein Embedded CA SC score</b>. It combines side-chain orientation and global distance based measures in an integrated framework using the united-residue model of polypeptide conformation for computing model-native similarity. SPECS captures both global and local quality aspects when evaluating structural similarity and is sensitive to minute variations in side-chain, thereby being a robust evaluation metric covering a wide range of modeling scenarios and various aspects of structural similarity. SPECS has the value in (0,1], where 1 indicates a perfect match between two structures.

## Requirement
1. Linux system: SPECS is tested on on x86_64 redhat linux system
2. GCC compiler

## Installation
You can run SPECS from <a href="http://watson.cse.eng.auburn.edu/SPECS/">SPECS web-server</a>. However, for large-scale benchmarking, we strongly recommend that you install and run SPECS locally.
SPECS is a stand-alone application. To install, you just need to download and compile as follows,
```
$ git clone https://github.com/Bhattacharya-Lab/SPECS.git
$ cd SPECS/src
$ g++ -o SPECS SPECS.cpp
```
## Usage
To run SPECS, type
```
$ ./SPECS
```
You will see the following output
```
**********************************************************************
*                            SPECS                                   *
*          Superposition-based Protein Embedded CA SC score          *
* Range of SPECS:                                                    *
*     0.0 <= SPECS <= 1.0, higher scores indicate better similarity  *
* For comments, please email to bhattacharyad@auburn.edu             *
**********************************************************************


Error! Model file must be provided

Usage: ./SPECS -m model -n native
   -m model  : model pdb file
   -n native : native pdb file
   -h help   : this message
```
<b>Example command to run SPECS</b>
```
$ ./SPECS -m ../example/sample_model.pdb -n ../example/sample_native.pdb
```
## Output
SPECS offers dynamic input validation and provideds instantaneous output. Upon running the above example command, you will see following output
```
**********************************************************************
*                           SPECS-SCORE                              *
*          Superposition-based Protein Embedded CA SC score          *
* Range of SPECS:                                                    *
*     0.0 <= SPECS <= 1.0, higher scores indicate better similarity  *
* For comments, please email to bhattacharyad@auburn.edu             *
**********************************************************************

Structure1: ../example/sample_model.pdb Length = 75
Structure2: ../example/sample_native.pdb Length = 62
Number of residues in common = 62
RMSD of common residues = 24.227

SPECS-score  = 0.3371 (dCA = 0.4395 rSC = 0.2921 angl = 0.2532 ang2 = 0.2790 tors = 0.1144)
TM-score     = 0.4688 (d0 = 2.67)
MaxSub-score = 0.4627 (d0 = 3.50)
GDT-TS-score = 0.5161 %(d<1) = 0.4355 %(d<2) = 0.4839 %(d<4) = 0.5484 %(d<8) = 0.5968
GDT-HA-score = 0.4395 %(d<0.5) = 0.2903 %(d<1) = 0.4355 %(d<2) = 0.4839 %(d<4) = 0.5484

-------- rotation matrix to rotate Chain-1 to Chain-2 ------
i             t(i)                u(i,1)               u(i,2)               u(i,3)
1        -3.3315000534        -0.2375999987         0.5454000235        -0.8037999868
2        14.3242998123         0.4174999893        -0.6898000240        -0.5914999843
3        64.4764022827        -0.8769999743        -0.4762000144        -0.0637999997

Superposition in the TM-score: Length(d<5.0) = 31 RMSD = 1.12
(':' denotes the residue pairs of distance < 5.0 Amstrong)
GHMEGKPKMEPAASSQAAVEELRTQVRELRSIIETMKDQQKREIKQLLSELDEEKKIRLRLQMEVNDIKKALQSK
              :::::::::::::::::::::::::::::::                              
----------PAASSQAAVEELRTQVRELRSIIETMKDQQKREIKQLLSELDEEKKIRLRLQMEVNDIKKAL---
123456789012345678901234567890123456789012345678901234567890123456789012345
```
## Contact
For questions and comments, please contact,<br/>
bhattacharyad@auburn.edu

## Reference
Alapati, R., Shuvo, MH., Bhattacharya, D. (2019) SPECS: Integration of side-chain orientation and global distance-based measures for improved evaluation of protein structural models. Submitted.
