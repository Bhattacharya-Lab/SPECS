# SPECS

SPECS is a side-chain-orientation-included protein model-native similarity metric for improved evaluation of protein structural models. SPECS stands for <b>Superposition-based Protein Embedded CA SC</b> score. It combines side-chain orientation and global distance based measures in an integrated framework using the united-residue model of polypeptide conformation for computing model-native similarity. SPECS captures both global and local quality aspects when evaluating structural similarity and is sensitive to minute variations in side-chain, thereby being a robust evaluation metric covering a wide range of modeling scenarios and various aspects of structural similarity. SPECS has the value in (0,1], where 1 indicates a perfect match between two structures.

If you find SPECS useful, please cite our <a href="https://doi.org/10.1371/journal.pone.0228245">PLOS ONE paper</a>.

## Requirement
1. Linux system: SPECS was tested on 64-bit Linux system
2. GCC compiler

## Installation
You can run SPECS from <a href="http://watson.cse.eng.auburn.edu/SPECS/">SPECS web-server</a>. However, for large-scale benchmarking, we strongly recommend that you install and run SPECS locally.
SPECS is a standalone application. To install, you just need to download and compile as follows:
```
$ git clone https://github.com/Bhattacharya-Lab/SPECS.git
$ cd SPECS/src
$ g++ -o SPECS SPECS.cpp
```
## Usage
To run SPECS, type
```
$ ./SPECS -h
```
You should see the following output:
```
**********************************************************************
*                            SPECS                                   *
*          Superposition-based Protein Embedded CA SC score          *
* Range of SPECS:                                                    *
*     0.0 <= SPECS <= 1.0, higher scores indicate better similarity  *
* For comments, please email to bhattacharyad@auburn.edu             *
**********************************************************************

Usage: ./SPECS -m model -n native
   -m model  : model pdb file
   -n native : native pdb file
   -h help   : this message
```
<b>Example command to run SPECS</b>
```
$ ./SPECS -m ../example/model.pdb -n ../example/native.pdb
```
## Output
SPECS provides instantaneous output and computes other structural similarity metrics. For ease of parsing, the output format of SPECS has been kept similar to the currently popular <a href="https://zhanglab.ccmb.med.umich.edu/TM-score/">TM-score</a> program. Upon running the above example command, you should see following output:
```
**********************************************************************
*                           SPECS-SCORE                              *
*          Superposition-based Protein Embedded CA SC score          *
* Range of SPECS:                                                    *
*     0.0 <= SPECS <= 1.0, higher scores indicate better similarity  *
* For comments, please email to bhattacharyad@auburn.edu             *
**********************************************************************

Structure1: ../example/model.pdb Length = 72
Structure2: ../example/native.pdb Length = 72
Number of residues in common = 72
RMSD of common residues = 3.323

SPECS-score  = 0.3704 (dCA = 0.4722 rSC = 0.3823 angl = 0.1597 ang2 = 0.3903 tors = 0.1417)
TM-score     = 0.6443 (d0 = 2.97)
MaxSub-score = 0.6197 (d0 = 3.50)
GDT-TS-score = 0.6840 %(d<1) = 0.3611 %(d<2) = 0.5556 %(d<4) = 0.8333 %(d<8) = 0.9861
GDT-HA-score = 0.4722 %(d<0.5) = 0.1389 %(d<1) = 0.3611 %(d<2) = 0.5556 %(d<4) = 0.8333

-------- rotation matrix to rotate Chain-1 to Chain-2 ------
i             t(i)                u(i,1)               u(i,2)               u(i,3)
1       -23.6888008118         0.7470999956         0.6341999769         0.1993000060
2        12.4891996384        -0.0847999975        -0.2064999938         0.9747999907
3        45.8860015869         0.6593000293        -0.7451000214        -0.1005000025

Superposition in the TM-score: Length(d<5.0) = 62 RMSD = 2.26
(':' denotes the residue pairs of distance < 5.0 Amstrong)
RLALSDAHFRRICQLIYQRAGIVLADHKRDMVYNRLVRRLRALGLDDFGRYLSMLEANQNSAEWQAFINALT
::::::::::::::::::::::    :::::::::::::::::::::::::::::::::     :::: :::
RLALSDAHFRRICQLIYQRAGIVLADHKRDMVYNRLVRRLRALGLDDFGRYLSMLEANQNSAEWQAFINALT
123456789012345678901234567890123456789012345678901234567890123456789012
```
