/*************************************************************************
 *                           SPECS
 * Superposition-based Protein Embedded CA SC score
 *
 * Developed by: Rahul Alapati
 * This program is used to compare two protein structures by applying
 * superposition-based integrated evaluation of protein structural models 
 * based on the distances of CA atoms as well as the distances and orientations 
 * of side-chain atoms to compute model versus native structural similarity.
 * Both structures must be in the PDB format. Users may obtain a brief
 * instruction by running the program in help mode, with -h argument.
 * For comments, please email to bhattacharyad@auburn.edu.
 * 
 * Copyright (C) 2019 Debswapna Bhattacharya. All rights reserved.
 * Permission to use, copy, modify, and distribute this program for
 * any purpose, with or without fee, is hereby granted, provided that
 * the notices on the head, the reference information, and this
 * copyright notice appear in all copies or substantial portions of
 * the Software. It is provided "as is" without express or implied
 * warranty.
 ************************** Change log ***********************************
 *     05/13/2019: The first version released.
 *************************************************************************/
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <complex>
#include <cmath>
#include <iomanip>

#define X 0
#define Y 1
#define Z 2

using namespace std;
using std::vector;

// amino acid residue types
const int ALA = 0;
const int CYS = 1;
const int ASP = 2;
const int GLU = 3;
const int PHE = 4;
const int GLY = 5;
const int HIS = 6;
const int ILE = 7;
const int LYS = 8;
const int LEU = 9;
const int MET = 10;
const int ASN = 11;
const int PRO = 12;
const int GLN = 13;
const int ARG = 14;
const int SER = 15;
const int THR = 16;
const int VAL = 17;
const int TRP = 18;
const int TYR = 19;

// amino acide residue code to three letter format
string seq3[20] = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"};

// amino acide residue code to one letter format
string seq[20] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};
 
//point3d object
struct point3d {
    long double x;
    long double y;
    long double z;
};

//pdbInfo object
struct pdbInfo {
    int id;
    int aa;
    char seq[3];
    point3d ca;
    point3d ca_05;
    point3d ca_1;
    point3d ca_2;
    point3d ca_4;
    point3d sc;
};

struct struc {
	long double xt[3001];
	long double yt[3001];
	long double zt[3001];
}s1; 

struct ca_struc {
        long double xt[3001];
        long double yt[3001];
        long double zt[3001];
}ca, ca_05, ca_1, ca_2, ca_4;

struct sc_struc {
        long double xt[3001];
        long double yt[3001];
        long double zt[3001];
}sc;

//align object
struct struct_align {
    int n_ali;
    int iA[3000];
    int iB[3000];
}align;

//global variables
char modFile[200] = "";
char natFile[200] = "";
float d0_fix = 0.0;
float ten_fix = 0.0;
char outname[200] = "";
string modSeq = " ";
string natSeq = " ";
int m_out = -1;
int m_fix = -1;
int m_len = -1;
float ratio = 1.0;
long double r_1[4][3000];
long double r_2[4][3000];
long double w[3000];
long double u[4][4];
long double t[4];
long double iq[3000];
long double rms;
long double armsd;
long double score;
float d, dca, dca_05, dca_1, dca_2, dca_4, rsc;
float d0;
int L_ini[3000];
int k_ali[3000];
int i_ali[3000];
int k_ali0[3000];
int ier = 0;
int ka0 = 0;
int n_cut = 0;
long double score_max;   //TMScore
long double score_maxsub_max; //MaxSub Score
long double score10_max; //TMScore10
double n_GDT05_max; //no of residues < 0.5
double n_GDT1_max; //no of residues < 1
double n_GDT2_max; //no of residues < 2
double n_GDT4_max; //no of residues < 4
double n_GDT8_max; //no of residues < 8
long double score_maxsub,score10;
double n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8;
string sequenceA = " ";
string sequenceB = " ";
string sequenceM = " ";
float d_output;

//for SPECS
long double U[4][4], U_SP_05[4][4], U_SP_1[4][4], U_SP_2[4][4], U_SP_4[4][4];
long double T[4], T_SP_05[4], T_SP_1[4], T_SP_2[4], T_SP_4[4];

//indicators
bool mFile = false;
bool nFile = false;

//function prototypes
void parseNextItem(int argc, char ** argv, int & i);
void parseCommandLine(int argc, char ** argv);
int getAA(const char * aa);
string loadPdb(char *filename, vector<pdbInfo> &pdb);
void calScores(vector<pdbInfo> &modPdb, vector<pdbInfo> &natPdb, long double &SPECS, float &d_GDT_HA, float &r_GDC_SC, float &theta_1, float &theta_2, float &phi, long double &TMScore, long double &RMSD, int &LCOMM, float &d0_print, long double &maxsubscore, int &GDT05, int &GDT1, int &GDT2, int &GDT4, int &GDT8, long double t[4], long double u[4][4], int &area_ctr, int &d_ctr);
long double u3b(long double w[3000], long double x[4][3000], long double y[4][3000], int n, int mode, long double rms, long double u[4][4], long double t[4], int ier);
void score_fun(struc s,vector<pdbInfo> &natPdb);
string NameMap(const char * residue);
point3d getDifference(point3d & p1, point3d &p2);
double getDotProduct(point3d & p1, point3d &p2);
point3d getUnit(point3d & p);
double getNorm(point3d & p);
point3d getCentroid(vector<point3d> &pointCloud);

/*************************************************************************
 * Name        : main
 * Purpose     : Compute SPECS-score
 * Arguments   : int argc, char ** argv
 * Return Type : int
 *************************************************************************/
int main(int argc, char ** argv) {

    //display header info
    cout << endl;
    cout << "**********************************************************************" << endl;
    cout << "*                           SPECS-SCORE                              *" << endl;
    cout << "*          Superposition-based Protein Embedded CA SC score          *" << endl; 
    cout << "* Range of SPECS:                                                    *" << endl;
    cout << "*     0.0 <= SPECS <= 1.0, higher scores indicate better similarity  *" << endl;
    cout << "* For comments, please email to bhattacharyad@auburn.edu             *" << endl;
    cout << "**********************************************************************" << endl;
    cout << endl;	
    //parse commandline and get the required arguements
    parseCommandLine(argc, argv);
    
    //load model
    vector<pdbInfo> modPdb;
    modSeq = loadPdb(modFile, modPdb);
    modSeq = modSeq.substr(1,modSeq.length());
   
    //load native
    vector<pdbInfo> natPdb;
    natSeq = loadPdb(natFile, natPdb);
    natSeq = natSeq.substr(1,natSeq.length());
 
    long double U[4][4];
    long double T[4];

    //calculate different similarity scores
    long double SPECS=1.0;
    float dCA=0.0, dSC=0.0, ang1=0.0, ang2=0.0, tors=0.0;
    long double TMSCORE=1.0;
    long double RMSD=1.0;
    int LCOMM=10;
    float d0_print = 0.0;
    long double maxsubscore = 0.0;
    int GDT05 = 0;
    int GDT1 = 0;
    int GDT2 = 0;
    int GDT4 = 0; 
    int GDT8 = 0;
    int area_ctr = 0;
    int d_ctr = 0;
    calScores(modPdb,natPdb,SPECS,dCA,dSC,ang1,ang2,tors,TMSCORE,RMSD,LCOMM,d0_print,maxsubscore,GDT05,GDT1,GDT2,GDT4,GDT8,T,U,area_ctr,d_ctr);
    
    if (LCOMM == 0) {
	cout << "There are no common residues in the input structures" << endl;
	return EXIT_SUCCESS;
    }
    else {
    	cout << "Structure1: " << modFile << " Length = " << modPdb.size() << endl;
    	cout << "Structure2: " << natFile << " Length = " << natPdb.size() << endl;
    	if (m_len == 1) {
		cout << "TM-score is normalized by " << ten_fix << endl;
    	}
    	cout << "Number of residues in common = " << LCOMM << endl;   
	cout << fixed << setprecision(3); 
    	cout << "RMSD of common residues = " << roundf(RMSD * 1000)/1000 << endl;
    	cout << endl;
	cout << fixed << setprecision(4);
    	if (m_len == 1) {
		TMSCORE = TMSCORE * float(natPdb.size())/float(ten_fix);
		SPECS = (SPECS/float(natPdb.size()-1)) * float(natPdb.size()-1)/float(ten_fix);	
    	}
  	
    	cout << "SPECS-score  = " << roundf(SPECS * 10000)/10000 << " (dCA = " << dCA << " rSC = " << dSC << " angl = " << ang1 << " ang2 = " << ang2 << " tors = " << tors << ")" << endl;
    	cout << "TM-score     = " << roundf(TMSCORE * 10000)/10000 << " (d0 = " << setprecision(2) << roundf(d0_print * 100)/100 << ")" << endl;
	cout << fixed << setprecision(4);
    	cout << "MaxSub-score = " << roundf(maxsubscore * ratio * 10000)/10000 << " (d0 = 3.50)" << endl;
    	float score_GDT = (GDT1 + GDT2 + GDT4 + GDT8)/float(4 * natPdb.size());
    	cout << "GDT-TS-score = " << roundf(score_GDT * ratio * 10000)/10000 << " %(d<1) = " << roundf(GDT1/float(natPdb.size()) * ratio * 10000)/10000 << " %(d<2) = " << roundf(GDT2/float(natPdb.size()) * ratio * 10000)/10000 << " %(d<4) = " << roundf(GDT4/float(natPdb.size()) * ratio * 10000)/10000 << " %(d<8) = " << roundf(GDT8/float(natPdb.size()) * ratio * 10000)/10000 << endl;
    	float score_GDT_HA = (GDT05 + GDT1 + GDT2 + GDT4)/float(4 * natPdb.size());
    	cout << "GDT-HA-score = " << roundf(score_GDT_HA * ratio * 10000)/10000 << " %(d<0.5) = " << roundf(GDT05/float(natPdb.size()) * ratio * 10000)/10000  << " %(d<1) = " << roundf(GDT1/float(natPdb.size()) * ratio * 10000)/10000 << " %(d<2) = " << roundf(GDT2/float(natPdb.size()) * ratio * 10000)/10000 << " %(d<4) = " << roundf(GDT4/float(natPdb.size()) * ratio * 10000)/10000 << endl;
    	cout << endl;
    	cout << "-------- rotation matrix to rotate Chain-1 to Chain-2 ------" << endl;
    	cout << "i             t(i)                u(i,1)               u(i,2)               u(i,3)" << endl;
    	for (int i = 1; i <= 3; i++) {
    		cout << i << "       "; 
		cout << fixed << setprecision(10) << setw(14) << setfill(' ') << roundf(T[i] * 10000)/10000 << "       "; 
		cout << fixed << setprecision(10) << setw(14) << setfill(' ') << roundf(U[i][1] * 10000)/10000 << "       "; 
		cout << fixed << setprecision(10) << setw(14) << setfill(' ') << roundf(U[i][2] * 10000)/10000 << "       "; 
		cout << fixed << setprecision(10) << setw(14) << setfill(' ') << roundf(U[i][3] * 10000)/10000 << endl;
    	}
	cout << endl;
	cout << fixed << setprecision(2);
    	cout << "Superposition in the TM-score: Length(d<" << setprecision(1) << d_output << ") = " << n_cut << " RMSD = " << setprecision(2) << armsd << endl;
    	cout << "(':' denotes the residue pairs of distance < " << setprecision(1) << d_output << " Amstrong)" << endl;
    	cout << sequenceA << endl;
    	cout << sequenceM << endl;
    	cout << sequenceB << endl;
    	for (int i = 1; i <= sequenceA.length(); i++) {
    		cout << i % 10;
    	}
    	cout << endl;
    	return EXIT_SUCCESS;
    }	
}

/*************************************************************************
 * Name        : parseNextItem
 * Purpose     : parse next item in command line argument
 * Arguments   : int argc, char ** argv, int & i
 * Return Type : void
 *************************************************************************/
void parseNextItem(int argc, char ** argv, int & i) {
    if (strncmp(argv[i], "-m", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No model file provided" << endl << endl;
            cout << "Usage: " << argv[0] << " -m model -n native" << endl;
            cout << "   -m model  : model pdb file" << endl;
            cout << "   -n native : native pdb file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        mFile = true;
        strcpy(modFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-n", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No native file provided" << endl << endl;
            cout << "Usage: " << argv[0] << " -m model -n native" << endl;
            cout << "   -m model  : model pdb file" << endl;
            cout << "   -n native : native pdb file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        nFile = true;
        strcpy(natFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-d", 2) == 0) {
        m_fix = 1;
        d0_fix = atof(argv[++i]);
    }
    else if (strncmp(argv[i], "-l", 2) == 0) {
        m_len = 1;
        ten_fix = atof(argv[++i]);
    }
    else if (strncmp(argv[i], "-o", 2) == 0) {
        m_out = 1;
        strcpy(outname, argv[++i]);
    }
    else if (strncmp(argv[i], "-h", 2) == 0) {
            cout << endl;
            cout << "Usage: " << argv[0] << " -m model -n native" << endl;
            cout << "   -m model  : model pdb file" << endl;
            cout << "   -n native : native pdb file" << endl;
	    cout << "   -h help   : this message" << endl;
            exit(0);
    }
    else {
            cout << endl;
            cout << "Error! Invalid option" << endl << endl;
            cout << "Usage: " << argv[0] << " -m model -n native" << endl;
            cout << "   -m model  : model pdb file" << endl;
            cout << "   -n native : native pdb file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);

    }
    i++;
}

/*************************************************************************
 * Name        : parseCommandLine
 * Purpose     : parse command line arguments
 * Arguments   : int argc, char ** argv
 * Return Type : void
 *************************************************************************/
void parseCommandLine(int argc, char ** argv) {
    int i = 1;
    while (i < argc)
        parseNextItem(argc, argv, i);
    if (!mFile) {
        cout << endl;
        cout << "Error! Model file must be provided" << endl << endl;
                    cout << "Usage: " << argv[0] << " -m model -n native" << endl;
                    cout << "   -m model  : model pdb file" << endl;
                    cout << "   -n native : native pdb file" << endl;
                    cout << "   -h help   : this message" << endl;
                    exit(0);
    }
    else {
        ifstream fin( modFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Model file not present" << endl << endl;
                                    cout << "Usage: " << argv[0] << " -m model -n native" << endl;
                                    cout << "   -m model  : model pdb file" << endl;
                                    cout << "   -n native : native pdb file" << endl;
                                    cout << "   -h help   : this message" << endl;
                                    exit(0);
        }
    }
    if (!nFile) {
        cout << endl;
        cout << "Error! Native file must be provided" << endl << endl;
                    cout << "Usage: " << argv[0] << " -m model -n native" << endl;
                    cout << "   -m model  : model pdb file" << endl;
                    cout << "   -n native : native pdb file" << endl;
                    cout << "   -h help   : this message" << endl;
                    exit(0);
    }
    else {
        ifstream fin( natFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Native file not present" << endl << endl;
                                    cout << "Usage: " << argv[0] << " -m model -n native" << endl;
                                    cout << "   -m model  : model pdb file" << endl;
                                    cout << "   -n native : native pdb file" << endl;
                                    cout << "   -h help   : this message" << endl;
                                    exit(0);
        }
    }
}

/*************************************************************************
 * Name        : NameMap
 * Purpose     : get NameMap
 * Arguments   : const char * aa
 * Return Type : string
 *************************************************************************/
string NameMap(const char * aa) {
	if (strlen(aa) == 3) {
        	if (strcmp(aa, "ALA") == 0)
            		return "A";
        	else if (strcmp(aa, "ARG") == 0)
            		return "R";
        	else if (strcmp(aa, "ASN") == 0)
            		return "N";
        	else if (strcmp(aa, "ASP") == 0)
            		return "D";
        	else if (strcmp(aa, "CYS") == 0)
            		return "C";
        	else if (strcmp(aa, "GLN") == 0)
            		return "Q";
        	else if (strcmp(aa, "GLU") == 0)
            		return "E";
        	else if (strcmp(aa, "GLY") == 0)
            		return "G";
        	else if (strcmp(aa, "HIS") == 0)
            		return "H";
        	else if (strcmp(aa, "ILE") == 0)
            		return "I";
        	else if (strcmp(aa, "LEU") == 0)
            		return "L";
        	else if (strcmp(aa, "LYS") == 0)
            		return "K";
        	else if (strcmp(aa, "MET") == 0)
            		return "M";
         	else if (strcmp(aa, "PHE") == 0)
            		return "F";
        	else if (strcmp(aa, "PRO") == 0)
            		return "P";
        	else if (strcmp(aa, "SER") == 0)
            		return "S";
        	else if (strcmp(aa, "THR") == 0)
            		return "T";
		else if (strcmp(aa, "TRP") == 0)
            		return "W";
        	else if (strcmp(aa, "TYR") == 0)
            		return "Y";
        	else if (strcmp(aa, "VAL") == 0)
            		return "V";
        	else {
            		cout << "Error! Invalid amino acid " << aa << endl;
            		exit(0);
        	}
    	}
	else if (strlen(aa) == 1) {
        	if (aa[0] == 'A')
            		return "ALA";
        	else if (aa[0] == 'C')
            		return "CYS";
        	else if (aa[0] == 'D')
            		return "ASP";
         	else if (aa[0] == 'E')
            		return "GLU";
        	else if (aa[0] == 'F')
            		return "PHE";
        	else if (aa[0] == 'G')
            		return "GLY";
        	else if (aa[0] == 'H')
            		return "HIS";
        	else if (aa[0] == 'I')
            		return "ILE";
        	else if (aa[0] == 'K')
            		return "LYS";
        	else if (aa[0] == 'L')
            		return "LEU";
        	else if (aa[0] == 'M')
            		return "MET";
        	else if (aa[0] == 'N')
            		return "ASN";
        	else if (aa[0] == 'P')
            		return "PRO";
        	else if (aa[0] == 'Q')
            		return "GLN";
        	else if (aa[0] == 'R')
            		return "ARG";
        	else if (aa[0] == 'S')
            		return "SER";
		else if (aa[0] == 'T')
            		return "THR";
        	else if (aa[0] == 'V')
            		return "VAL";
        	else if (aa[0] == 'W')
            		return "TRP";
        	else if (aa[0] == 'Y')
            		return "TYR";
        	else {
            		cout << "Error! Invalid amino acid " << aa << endl;
            		exit(0);
        	}
    	}
    	else {
        	cout << "Error! Invalid amino acid " << aa << endl;
		exit(0);
	}
}

/*************************************************************************
 * Name        : getAA
 * Purpose     : convert AA name to a numerical code
 * Arguments   : const char * aa
 * Return Type : int
 *************************************************************************/
int getAA(const char * aa) {
    if (strlen(aa) == 3) {
        if (strcmp(aa, "ALA") == 0)
            return (ALA);
        else if (strcmp(aa, "ARG") == 0)
            return (ARG);
        else if (strcmp(aa, "ASN") == 0)
            return (ASN);
        else if (strcmp(aa, "ASP") == 0)
            return (ASP);
        else if (strcmp(aa, "CYS") == 0)
            return (CYS);
        else if (strcmp(aa, "GLN") == 0)
            return (GLN);
        else if (strcmp(aa, "GLU") == 0)
            return (GLU);
        else if (strcmp(aa, "GLY") == 0)
            return (GLY);
        else if (strcmp(aa, "HIS") == 0)
            return (HIS);
        else if (strcmp(aa, "ILE") == 0)
            return (ILE);
        else if (strcmp(aa, "LEU") == 0)
            return (LEU);
        else if (strcmp(aa, "LYS") == 0)
            return (LYS);
        else if (strcmp(aa, "MET") == 0)
            return (MET);
	 else if (strcmp(aa, "PHE") == 0)
            return (PHE);
        else if (strcmp(aa, "PRO") == 0)
            return (PRO);
        else if (strcmp(aa, "SER") == 0)
            return (SER);
        else if (strcmp(aa, "THR") == 0)
            return (THR);
        else if (strcmp(aa, "TRP") == 0)
            return (TRP);
        else if (strcmp(aa, "TYR") == 0)
            return (TYR);
        else if (strcmp(aa, "VAL") == 0)
            return (VAL);
        else {
            cout << "Error! Invalid amino acid " << aa << endl;
            exit(0);
        }
    }
    else if (strlen(aa) == 1) {
        if (aa[0] == 'A')
            return ALA;
        else if (aa[0] == 'C')
            return CYS;
        else if (aa[0] == 'D')
            return ASP;
	 else if (aa[0] == 'E')
            return GLU;
        else if (aa[0] == 'F')
            return PHE;
        else if (aa[0] == 'G')
            return GLY;
        else if (aa[0] == 'H')
            return HIS;
        else if (aa[0] == 'I')
            return ILE;
        else if (aa[0] == 'K')
            return LYS;
        else if (aa[0] == 'L')
            return LEU;
        else if (aa[0] == 'M')
            return MET;
        else if (aa[0] == 'N')
            return ASN;
        else if (aa[0] == 'P')
            return PRO;
        else if (aa[0] == 'Q')
            return GLN;
        else if (aa[0] == 'R')
            return ARG;
        else if (aa[0] == 'S')
            return SER;
        else if (aa[0] == 'T')
	    return THR;
        else if (aa[0] == 'V')
            return VAL;
        else if (aa[0] == 'W')
            return TRP;
        else if (aa[0] == 'Y')
            return TYR;
        else {
            cout << "Error! Invalid amino acid " << aa << endl;
            exit(0);
        }
    }
    else {
        cout << "Error! Invalid amino acid " << aa << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : loadPdb
 * Purpose     : loads a pdb file into a vector of pdbInfo object
 * Arguments   : char *filename, vector<pdbInfo> &pdb
 * Return Type : string
 *************************************************************************/
string loadPdb(char *filename, vector<pdbInfo> &pdb) {
    string line, str;
    string returnseq = " ";
    string atom ("ATOM ");
    int prevRes = -999999;
    point3d caAtom;
    point3d scAtom;
    vector<point3d> scAtomCloud;
    pdbInfo pdbData;
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            if(line.compare(0, atom.length(), atom)==0) {
		// consider the first alternate location id (i.e. A) if present
		if( line.compare(16, 1, " ") == 0 || line.compare(16, 1, "A") == 0 ) {
		    // get the CA atom coordinate
		    if( line.compare(12, 4, "CA  ") == 0 || line.compare(12, 4, " CA ") == 0 || line.compare(12, 4, "  CA") == 0 ) {	
                	int res = atoi(line.substr(22, 4).c_str());
                	int anmae = getAA(line.substr(17, 3).c_str());
			char seqA[3];
			strcpy(seqA,line.substr(17, 3).c_str());
			// seek for the next residue
			if (res != prevRes) {
                    		if (prevRes != -999999) {
					// insert the data collected so far
					pdbData.ca = caAtom;
					pdbData.sc = getCentroid(scAtomCloud);
                        		pdb.push_back(pdbData);
                    		}
                    		prevRes = res;
		    		scAtomCloud.clear();
                    		pdbData.id = res;
                    		pdbData.aa = anmae;
		    		strcpy(pdbData.seq,seqA);
                	}
                        returnseq = returnseq + NameMap(line.substr(17,3).c_str());
			caAtom.x = atof(line.substr(30, 8).c_str());
                        caAtom.y = atof(line.substr(38, 8).c_str());
                        caAtom.z = atof(line.substr(46, 8).c_str());
                    }
		    // get the sidechain heavy atoms
		    if(line.compare(12, 4, "N   ") != 0 && line.compare(12, 4, " N  ") != 0 && line.compare(12, 4, "  N ") != 0 &&
                       line.compare(12, 4, "C   ") != 0 && line.compare(12, 4, " C  ") != 0 && line.compare(12, 4, "  C ") != 0 &&
                       line.compare(12, 4, "O   ") != 0 && line.compare(12, 4, " O  ") != 0 && line.compare(12, 4, "  O ") != 0 &&
                       line.compare(12, 4, "CA  ") != 0 && line.compare(12, 4, " CA ") != 0 && line.compare(12, 4, "  CA") != 0) {
                        scAtom.x = atof(line.substr(30, 8).c_str());
                        scAtom.y = atof(line.substr(38, 8).c_str());
                        scAtom.z = atof(line.substr(46, 8).c_str());
                        scAtomCloud.push_back(scAtom);
                    }
                }
            }
        }
	fin.close();
	// for the last residue
	pdbData.ca = caAtom;
	pdbData.sc = getCentroid(scAtomCloud);
        pdb.push_back(pdbData);
	return returnseq;
    }
    else {
        cout << "Error! pdb file can not open " << filename << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : calScores
 * Purpose     : calculate different similarity scores
 * Arguments   : vector<pdbInfo> &modPdb, vector<pdbInfo> &natPdb
 * Return Type : long double &SPECS, long double &tmscore, long double &rmsd, int &lcomm, float &d0_print, long double &maxsubscore, int &GDT05, int &GDT1, int &GDT2, int &GDT4, int &GDT8,  *               long double T[4], long double U[4][4], int &area_ctr, int &d_ctr
 *************************************************************************/
void calScores(vector<pdbInfo> &modPdb, vector<pdbInfo> &natPdb, long double &SPECS, float &d_GDT_HA, float &r_GDC_SC, float &theta_1, float &theta_2, float &phi, long double &tmscore, long double &rmsd, int &lcomm, float &d0_print, long double &maxsubscore, int &GDT05, int &GDT1, int &GDT2, int &GDT4, int &GDT8, long double T[4], long double U[4][4], int &area_ctr, int &d_ctr)
{
	//initialize w with 1.0
	for (int i = 0; i < 3000; i++) {
		w[i] = 1.0;
	}

	//pickup aligned residues
	int k = 0;
	for (int i = 0; i < modPdb.size(); i++) {
		for (int j = 0; j < natPdb.size(); j++) {
			if (modPdb[i].id == natPdb[j].id) {
				k = k + 1;
				align.iA[k] = i;
				align.iB[k] = j;
				break;
			}
		}
	}
	
	//if no aligned residues
	align.n_ali = k;
	lcomm = align.n_ali;
	if (align.n_ali < 1) {
		tmscore = 0.0;
		rmsd = 0.0;
		return;
	}

	//compute d0 based on different conditions
	if (natPdb.size() > 15) {
		d0 = 1.24 * pow((natPdb.size() - 15),(1.0/3.0)) - 1.8;
	}
	else {
		d0 = 0.5;
	}
	
	if (m_len == 1) {
		d0 = 1.24 * pow((ten_fix - 15),(1.0/3.0)) - 1.8;		
	}
	
	if (m_fix == 1) {
		d0 = d0_fix;
	}

	if (d0 < 0.5) {
		d0 = 0.5;
	}
	
	//search for the optimal d0
	float d0_search = d0;
	if (d0_search > 8.0) {
		d0_search = 8.0;
	}

	if (d0_search < 4.5) {
		d0_search = 4.5;
	}

	//iterative parameters for searching the optimal d0, alignment and scores
	int n_it = 20; //max iterations
	d_output = 5; //for output alignment
	
	if (m_fix == 1) {
		d_output = d0_fix;
	}

	int n_init_max = 6; //max no of L_init
	int n_init = 0;
	int L_ini_min = 4;
	
	if (align.n_ali < 4) {
		L_ini_min = align.n_ali;
	}
	
	bool flag1 = false;
	for (int i = 1; i <= n_init_max-1; i++) {
		n_init = n_init + 1;
		L_ini[n_init] = (align.n_ali/(int) pow(2,n_init-1));
		
		if (L_ini[n_init] <= L_ini_min) {
			L_ini[n_init] = L_ini_min;
			flag1 = true;
			break;
		}
	}
	
	if (!flag1) {
		n_init = n_init + 1;
		L_ini[n_init] = L_ini_min;
	} 
	
	//iterative search strategy for finding the optimal d0, rotation and translation matrices, alignments and scores
	score_max = -1.0;   //TMScore
	score_maxsub_max = -1.0; //MaxSub Score
	score10_max = -1.0; //TMScore10
	n_GDT05_max = -1.0; //no of residues < 0.5
	n_GDT1_max = -1.0; //no of residues < 1
	n_GDT2_max = -1.0; //no of residues < 2
	n_GDT4_max = -1.0; //no of residues < 4	
	n_GDT8_max = -1.0; //no of residues < 8
	for (int i_init = 1; i_init <= n_init; i_init++) {
		int L_init = L_ini[i_init];
		int iL_max = align.n_ali - L_init + 1;
		for (int iL = 1; iL <= iL_max; iL++) {
			int LL = 0;
			int ka = 0;
			for (int i = 1; i <= L_init; i++) {
				k = iL + i - 1;
				r_1[1][i] = modPdb[align.iA[k]].ca.x;
				r_1[2][i] = modPdb[align.iA[k]].ca.y;
				r_1[3][i] = modPdb[align.iA[k]].ca.z;
				r_2[1][i] = natPdb[align.iB[k]].ca.x;
				r_2[2][i] = natPdb[align.iB[k]].ca.y;
				r_2[3][i] = natPdb[align.iB[k]].ca.z;
				ka = ka + 1;
				k_ali[ka] = k;
				LL = LL + 1;
			}

			if (i_init == 1) {
				//calculate rotation and translation matrices for an alignment
				rms = u3b(w,r_1,r_2,LL,2,rms,u,t,ier);
				armsd = sqrt(rms/LL);
				rmsd = armsd;
			}
			else {
				rms = u3b(w,r_1,r_2,LL,1,rms,u,t,ier);
			}
			
			for (int j = 0; j <= modPdb.size()-1; j++) {
				s1.xt[j] = (t[1] + u[1][1] * modPdb[j].ca.x + u[1][2] * modPdb[j].ca.y + u[1][3] * modPdb[j].ca.z);
				s1.yt[j] = (t[2] + u[2][1] * modPdb[j].ca.x + u[2][2] * modPdb[j].ca.y + u[2][3] * modPdb[j].ca.z);
				s1.zt[j] = (t[3] + u[3][1] * modPdb[j].ca.x + u[3][2] * modPdb[j].ca.y + u[3][3] * modPdb[j].ca.z);
			}
			
			d = d0_search - 1;
			//calculate scores based on the current alignment
			score_fun(s1, natPdb);
		        if (score_max < score) {
				score_max = score;
				ka0 = ka;
				for (int i = 1; i <= ka0; i++) {
					k_ali0[i] = k_ali[i];
				} 
			}

			if (score10_max < score10) {
				score10_max = score10;
			}
			if (score_maxsub_max < score_maxsub) {
				score_maxsub_max = score_maxsub;
			}
			if (n_GDT05_max < n_GDT05) {
				n_GDT05_max =  n_GDT05;

				//getting the optimal rotation and translation matrices for GDTHA threshold 0.5 A
				for (int i = 0; i < 4; i++) {
                                        T_SP_05[i] = t[i];
                                }

                                for (int i = 0; i < 4; i++) {
                                        for (int j = 0; j < 4; j++) {
                                                U_SP_05[i][j] = u[i][j];
                                        }
                                }
			}
			if (n_GDT1_max < n_GDT1) {
                                n_GDT1_max =  n_GDT1;

				//getting the optimal rotation and translation matrices for GDTHA threshold 1 A
				for (int i = 0; i < 4; i++) {
                                        T_SP_1[i] = t[i];
                                }

                                for (int i = 0; i < 4; i++) {
                                        for (int j = 0; j < 4; j++) {
                                                U_SP_1[i][j] = u[i][j];
                                        }
                                }
                        }
			if (n_GDT2_max < n_GDT2) {
                                n_GDT2_max =  n_GDT2;

				//getting the optimal rotation and translation matrices for GDTHA threshold 2 A
				for (int i = 0; i < 4; i++) {
                                        T_SP_2[i] = t[i];
                                }

                                for (int i = 0; i < 4; i++) {
                                        for (int j = 0; j < 4; j++) {
                                                U_SP_2[i][j] = u[i][j];
                                        }
                                }
                        }
			if (n_GDT4_max < n_GDT4) {
                                n_GDT4_max =  n_GDT4;
	
				//getting the optimal rotation and translation matrices for GDTHA threshold 4 A
				for (int i = 0; i < 4; i++) {
                                        T_SP_4[i] = t[i];
                                }

                                for (int i = 0; i < 4; i++) {
                                        for (int j = 0; j < 4; j++) {
                                                U_SP_4[i][j] = u[i][j];
                                        }
                                }
                        }
			if (n_GDT8_max < n_GDT8) {
                                n_GDT8_max =  n_GDT8;
                        }
			
			//final optimal search for computing the final scores
			d = d0_search + 1;
			for (int it = 1; it <= n_it; it++) {
				LL = 0;
				ka = 0;
				for (int i = 1; i <= n_cut; i++) {
					int m = i_ali[i];
					r_1[1][i] = modPdb[align.iA[m]].ca.x;
					r_1[2][i] = modPdb[align.iA[m]].ca.y;
					r_1[3][i] = modPdb[align.iA[m]].ca.z;
					r_2[1][i] = natPdb[align.iB[m]].ca.x;
					r_2[2][i] = natPdb[align.iB[m]].ca.y;
					r_2[3][i] = natPdb[align.iB[m]].ca.z;
					ka++;
					k_ali[ka] = m;
					LL++;		
				}	
				
				//calculate rotation and translation matrices for an alignment
				rms = u3b(w,r_1,r_2,LL,1,rms,u,t,ier);

				for (int j = 0; j <= modPdb.size()-1; j++) {
					s1.xt[j] = (t[1] + u[1][1] * modPdb[j].ca.x + u[1][2] * modPdb[j].ca.y + u[1][3] * modPdb[j].ca.z);
					s1.yt[j] = (t[2] + u[2][1] * modPdb[j].ca.x + u[2][2] * modPdb[j].ca.y + u[2][3] * modPdb[j].ca.z);
					s1.zt[j] = (t[3] + u[3][1] * modPdb[j].ca.x + u[3][2] * modPdb[j].ca.y + u[3][3] * modPdb[j].ca.z);
				}
				
				//calculate scores based on the current alignment
				score_fun(s1, natPdb);
				if (isnan(score)) 
					score = 0.0;
				if (score_max < score) {
					score_max = score;
					ka0 = ka;
					for (int i = 1; i <= ka; i++) {
						k_ali0[i] = k_ali[i];
					}
				}
				
				if (score10_max < score10) {
                                	score10_max = score10;
                        	}

                        	if (score_maxsub_max < score_maxsub) {
                                	score_maxsub_max = score_maxsub;
                        	}

                        	if (n_GDT05_max < n_GDT05) {
                                	n_GDT05_max =  n_GDT05;

					//getting the optimal rotation and translation matrices for GDTHA threshold 0.5 A
					for (int i = 0; i < 4; i++) {
                                        	T_SP_05[i] = t[i];
                                	}

                                	for (int i = 0; i < 4; i++) {
                                        	for (int j = 0; j < 4; j++) {
                                                	U_SP_05[i][j] = u[i][j];
                                        	}
                                	}
                        	}

                        	if (n_GDT1_max < n_GDT1) {
                                	n_GDT1_max =  n_GDT1;
	
					//getting the optimal rotation and translation matrices for GDTHA threshold 1 A
					for (int i = 0; i < 4; i++) {
                                                T_SP_1[i] = t[i];
                                        }

                                        for (int i = 0; i < 4; i++) {
                                                for (int j = 0; j < 4; j++) {
                                                        U_SP_1[i][j] = u[i][j];
                                                }
                                        }
                        	}

                        	if (n_GDT2_max < n_GDT2) {
                                	n_GDT2_max =  n_GDT2;

					//getting the optimal rotation and translation matrices for GDTHA threshold 2 A
					for (int i = 0; i < 4; i++) {
                                                T_SP_2[i] = t[i];
                                        }

                                        for (int i = 0; i < 4; i++) {
                                                for (int j = 0; j < 4; j++) {
                                                        U_SP_2[i][j] = u[i][j];
                                                }
                                        }
                        	}

                        	if (n_GDT4_max < n_GDT4) {
                                	n_GDT4_max =  n_GDT4;

					//getting the optimal rotation and translation matrices for GDTHA threshold 4 A
					for (int i = 0; i < 4; i++) {
                                                T_SP_4[i] = t[i];
                                        }

                                        for (int i = 0; i < 4; i++) {
                                                for (int j = 0; j < 4; j++) {
                                                        U_SP_4[i][j] = u[i][j];
                                                }
                                        }
                        	}

                        	if (n_GDT8_max < n_GDT8) {
                                	n_GDT8_max =  n_GDT8;
                        	}
				
				if (it == n_it) {
					break;
				}
				
				if (n_cut == ka) {
					int neq = 0;
					for (int i = 1; i <= n_cut; i++) {
						if (i_ali[i] == k_ali[i]) {
							neq = neq + 1;
						}
					}
					
					if (n_cut == neq) {
						break;
					}
				}
			}		
		}
	}
	
	ratio = 1;
	if (m_len > 0) {
		ratio = float(natPdb.size())/float(ten_fix);
	}
	
	int LL = 0;
	for (int i = 1; i <= ka0; i++) {
		int m = k_ali0[i];
		r_1[1][i] = modPdb[align.iA[m]].ca.x;
		r_1[2][i] = modPdb[align.iA[m]].ca.y;
		r_1[3][i] = modPdb[align.iA[m]].ca.z;
		r_2[1][i] = natPdb[align.iB[m]].ca.x;
		r_2[2][i] = natPdb[align.iB[m]].ca.y;
		r_2[3][i] = natPdb[align.iB[m]].ca.z;
		LL = LL + 1;
	}	 
	
	//calculate rotation and translation matrices for an alignment
	rms = u3b(w,r_1,r_2,LL,1,rms,u,t,ier);
	
	for (int j = 0; j <= modPdb.size()-1; j++) {
		s1.xt[j] = (t[1] + u[1][1] * modPdb[j].ca.x + u[1][2] * modPdb[j].ca.y + u[1][3] * modPdb[j].ca.z);
                s1.yt[j] = (t[2] + u[2][1] * modPdb[j].ca.x + u[2][2] * modPdb[j].ca.y + u[2][3] * modPdb[j].ca.z);
                s1.zt[j] = (t[3] + u[3][1] * modPdb[j].ca.x + u[3][2] * modPdb[j].ca.y + u[3][3] * modPdb[j].ca.z);
	}	
	
	//final d0, TMScore, MaxSub Score and the GDTHA thresholds
	d0_print = d0;
	tmscore = score_max;
	maxsubscore = score_maxsub_max;
        GDT05 = n_GDT05_max;
	GDT1 = n_GDT1_max;
	GDT2 = n_GDT2_max;
	GDT4 = n_GDT4_max;
	GDT8 = n_GDT8_max;

	//calculating rmsd in only superposed regions
	d = d_output;
	score_fun(s1,natPdb);
	
	for (int i = 0; i < modPdb.size(); i++) {
		iq[i] = 0.0;
	}
	
	//calculate distances based on aligned residues
	long double dis = 0.0;
	for (int i = 1; i <= n_cut; i++) {
		int j = align.iA[i_ali[i]];
		int k = align.iB[i_ali[i]];
		dis = sqrt((s1.xt[j] - natPdb[k].ca.x)*(s1.xt[j] - natPdb[k].ca.x) + (s1.yt[j] - natPdb[k].ca.y)*(s1.yt[j] - natPdb[k].ca.y) + (s1.zt[j] - natPdb[k].ca.z)*(s1.zt[j] - natPdb[k].ca.z));
		if (dis < d_output) {
			iq[j] = 1.0;
		}		
	}
	
	//for obtaining the sequence alignment
	k = 0;
	int i = 0;
	int j = 0;
	sequenceA = " ";
	sequenceB = " ";
	sequenceM = " ";
	
	while ((i <= modPdb.size()-1) || (j <= natPdb.size()-1)) {
		if ((i > modPdb.size()-1) && (j <= natPdb.size()-1)) {
			sequenceA = sequenceA + '-';
			sequenceB = sequenceB + natSeq[j];
			sequenceM = sequenceM + ' ';
			j++; 
		}
		else if ((i <= modPdb.size()-1) && (j > natPdb.size()-1)) {
			sequenceA = sequenceA + modSeq[i];
			sequenceB = sequenceB + '-';
			sequenceM = sequenceM + ' ';
                        i++;	
		}
		else if (modPdb[i].id == natPdb[j].id) {
			sequenceA = sequenceA + modSeq[i];
			sequenceB = sequenceB + natSeq[j];
			if (iq[i] == 1.0) {
				sequenceM = sequenceM + ':';
			}
			else {
				sequenceM = sequenceM + ' ';
			}
			i++;
			j++;		
		}
		else if (modPdb[i].id < natPdb[j].id) {
			if (i < modPdb.size()-1) {
				sequenceA = sequenceA + modSeq[i];
			}
			else {
				sequenceA = sequenceA + ' ';
			}
			sequenceB = sequenceB + '-';
			sequenceM = sequenceM + ' ';
			i++;
		}
		else if (modPdb[i].id > natPdb[j].id) {
			sequenceA = sequenceA + '-';
			if (j < natPdb.size()-1) {
				sequenceB = sequenceB + natSeq[j];
			}
			else {
				sequenceB = sequenceB + ' ';
			}
			sequenceM = sequenceM + ' ';
			j++;
		}
	}
	
	sequenceA = sequenceA.substr(1,sequenceA.size());
	sequenceB = sequenceB.substr(1,sequenceB.size());
	sequenceM = sequenceM.substr(1,sequenceM.size());	
	
	long double r_3[4][3000];
	long double r_4[4][3000];
	
	LL = 0;
        for (int i = 1; i <= n_cut; i++) {
         	int m = i_ali[i];
                r_3[1][i] = modPdb[align.iA[m]].ca.x;
                r_3[2][i] = modPdb[align.iA[m]].ca.y;
                r_3[3][i] = modPdb[align.iA[m]].ca.z;
                r_4[1][i] = natPdb[align.iB[m]].ca.x;
                r_4[2][i] = natPdb[align.iB[m]].ca.y;
                r_4[3][i] = natPdb[align.iB[m]].ca.z;
                LL++;
        }
	for (int i = 0; i < 4; i++) {
		T[i] = t[i];
	}
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			U[i][j] = u[i][j];
		}
	}
	
	rms = u3b(w,r_3,r_4,LL,0,rms,u,t,ier);
        armsd = sqrt(rms/LL);
	
	//translating and rotating ca and sc atoms based on the final t and u matrices	
	vector<pdbInfo> transPdb;
	pdbInfo pdbdata;
	point3d caatom, caatom_05, caatom_1, caatom_2, caatom_4;
    	point3d scatom;
	for (int j = 0; j <= modPdb.size()-1; j++) {
		//for SPECS CA calculation
                ca.xt[j] = (T[1] + U[1][1] * modPdb[j].ca.x + U[1][2] * modPdb[j].ca.y + U[1][3] * modPdb[j].ca.z);
                ca.yt[j] = (T[2] + U[2][1] * modPdb[j].ca.x + U[2][2] * modPdb[j].ca.y + U[2][3] * modPdb[j].ca.z);
                ca.zt[j] = (T[3] + U[3][1] * modPdb[j].ca.x + U[3][2] * modPdb[j].ca.y + U[3][3] * modPdb[j].ca.z);
		
		//for GDTHA 0.5 A threshold calculation
		ca_05.xt[j] = (T_SP_05[1] + U_SP_05[1][1] * modPdb[j].ca.x + U_SP_05[1][2] * modPdb[j].ca.y + U_SP_05[1][3] * modPdb[j].ca.z);
                ca_05.yt[j] = (T_SP_05[2] + U_SP_05[2][1] * modPdb[j].ca.x + U_SP_05[2][2] * modPdb[j].ca.y + U_SP_05[2][3] * modPdb[j].ca.z);
                ca_05.zt[j] = (T_SP_05[3] + U_SP_05[3][1] * modPdb[j].ca.x + U_SP_05[3][2] * modPdb[j].ca.y + U_SP_05[3][3] * modPdb[j].ca.z);

		//for GDTHA 1 A threshold calculation
		ca_1.xt[j] = (T_SP_1[1] + U_SP_1[1][1] * modPdb[j].ca.x + U_SP_1[1][2] * modPdb[j].ca.y + U_SP_1[1][3] * modPdb[j].ca.z);
                ca_1.yt[j] = (T_SP_1[2] + U_SP_1[2][1] * modPdb[j].ca.x + U_SP_1[2][2] * modPdb[j].ca.y + U_SP_1[2][3] * modPdb[j].ca.z);
                ca_1.zt[j] = (T_SP_1[3] + U_SP_1[3][1] * modPdb[j].ca.x + U_SP_1[3][2] * modPdb[j].ca.y + U_SP_1[3][3] * modPdb[j].ca.z);

		//for GDTHA 2 A threshold calculation
		ca_2.xt[j] = (T_SP_2[1] + U_SP_2[1][1] * modPdb[j].ca.x + U_SP_2[1][2] * modPdb[j].ca.y + U_SP_2[1][3] * modPdb[j].ca.z);
                ca_2.yt[j] = (T_SP_2[2] + U_SP_2[2][1] * modPdb[j].ca.x + U_SP_2[2][2] * modPdb[j].ca.y + U_SP_2[2][3] * modPdb[j].ca.z);
                ca_2.zt[j] = (T_SP_2[3] + U_SP_2[3][1] * modPdb[j].ca.x + U_SP_2[3][2] * modPdb[j].ca.y + U_SP_2[3][3] * modPdb[j].ca.z);

		//for GDTHA 4 A threshold calculation
		ca_4.xt[j] = (T_SP_4[1] + U_SP_4[1][1] * modPdb[j].ca.x + U_SP_4[1][2] * modPdb[j].ca.y + U_SP_4[1][3] * modPdb[j].ca.z);
                ca_4.yt[j] = (T_SP_4[2] + U_SP_4[2][1] * modPdb[j].ca.x + U_SP_4[2][2] * modPdb[j].ca.y + U_SP_4[2][3] * modPdb[j].ca.z);
                ca_4.zt[j] = (T_SP_4[3] + U_SP_4[3][1] * modPdb[j].ca.x + U_SP_4[3][2] * modPdb[j].ca.y + U_SP_4[3][3] * modPdb[j].ca.z);
		
		//for SPECS SC calculation
		sc.xt[j] = (T[1] + U[1][1] * modPdb[j].sc.x + U[1][2] * modPdb[j].sc.y + U[1][3] * modPdb[j].sc.z);
                sc.yt[j] = (T[2] + U[2][1] * modPdb[j].sc.x + U[2][2] * modPdb[j].sc.y + U[2][3] * modPdb[j].sc.z);
                sc.zt[j] = (T[3] + U[3][1] * modPdb[j].sc.x + U[3][2] * modPdb[j].sc.y + U[3][3] * modPdb[j].sc.z);
		
		caatom.x = ca.xt[j];
		caatom.y = ca.yt[j];
		caatom.z = ca.zt[j];
		
		caatom_05.x = ca_05.xt[j];
                caatom_05.y = ca_05.yt[j];
                caatom_05.z = ca_05.zt[j];

		caatom_1.x = ca_1.xt[j];
                caatom_1.y = ca_1.yt[j];
                caatom_1.z = ca_1.zt[j];

		caatom_2.x = ca_2.xt[j];
                caatom_2.y = ca_2.yt[j];
                caatom_2.z = ca_2.zt[j];

		caatom_4.x = ca_4.xt[j];
                caatom_4.y = ca_4.yt[j];
                caatom_4.z = ca_4.zt[j];
	
		scatom.x = sc.xt[j];
		scatom.y = sc.yt[j];
		scatom.z = sc.zt[j];

		pdbdata.ca = caatom;
		pdbdata.ca_05 = caatom_05;
		pdbdata.ca_1 = caatom_1;
		pdbdata.ca_2 = caatom_2;
		pdbdata.ca_4 = caatom_4;
        	pdbdata.sc = scatom;
        	transPdb.push_back(pdbdata);
        }
	
	//final SPECS score calculation	
	double d_GDT05 = 0.0,d_GDT1 = 0.0,d_GDT2 = 0.0,d_GDT4 = 0.0,d_GDT8 = 0.0;
	double n_thetaij1_120 = 0.0,n_thetaij1_90 = 0.0,n_thetaij1_60 = 0.0,n_thetaij1_30 = 0.0;
	double n_thetaij2_120 = 0.0,n_thetaij2_90 = 0.0,n_thetaij2_60 = 0.0,n_thetaij2_30 = 0.0; 
	double n_phiij_300 = 0.0,n_phiij_270 = 0.0,n_phiij_240 = 0.0,n_phiij_210 = 0.0,n_phiij_180 = 0.0,n_phiij_150 = 0.0,n_phiij_120 = 0.0,n_phiij_90 = 0.0,n_phiij_60 = 0.0,n_phiij_30 = 0.0;
	double r_GDT5 = 0.0,r_GDT45 = 0.0,r_GDT4 = 0.0,r_GDT35 = 0.0,r_GDT3 = 0.0,r_GDT25 = 0.0,r_GDT2 = 0.0,r_GDT15 = 0.0,r_GDT1 = 0.0,r_GDT05 = 0.0;	
	for (int k = 1; k <= align.n_ali; k++) {
                        int i = align.iA[k]; // ![1,nseqA] reorder number of structure A
                        int j = align.iB[k]; // ![1,nseqB]
			//GDT-HA calculation using CA atoms
                        dca = sqrt((ca.xt[i] - natPdb[j].ca.x)*(ca.xt[i] - natPdb[j].ca.x) + (ca.yt[i] - natPdb[j].ca.y)*(ca.yt[i] - natPdb[j].ca.y) + (ca.zt[i] - natPdb[j].ca.z)*(ca.zt[i] - natPdb[j].ca.z));
			dca_05 = sqrt((ca_05.xt[i] - natPdb[j].ca.x)*(ca_05.xt[i] - natPdb[j].ca.x) + (ca_05.yt[i] - natPdb[j].ca.y)*(ca_05.yt[i] - natPdb[j].ca.y) + (ca_05.zt[i] - natPdb[j].ca.z)*(ca_05.zt[i] - natPdb[j].ca.z));
			dca_1 = sqrt((ca_1.xt[i] - natPdb[j].ca.x)*(ca_1.xt[i] - natPdb[j].ca.x) + (ca_1.yt[i] - natPdb[j].ca.y)*(ca_1.yt[i] - natPdb[j].ca.y) + (ca_1.zt[i] - natPdb[j].ca.z)*(ca_1.zt[i] - natPdb[j].ca.z));
			dca_2 = sqrt((ca_2.xt[i] - natPdb[j].ca.x)*(ca_2.xt[i] - natPdb[j].ca.x) + (ca_2.yt[i] - natPdb[j].ca.y)*(ca_2.yt[i] - natPdb[j].ca.y) + (ca_2.zt[i] - natPdb[j].ca.z)*(ca_2.zt[i] - natPdb[j].ca.z));
			dca_4 = sqrt((ca_4.xt[i] - natPdb[j].ca.x)*(ca_4.xt[i] - natPdb[j].ca.x) + (ca_4.yt[i] - natPdb[j].ca.y)*(ca_4.yt[i] - natPdb[j].ca.y) + (ca_4.zt[i] - natPdb[j].ca.z)*(ca_4.zt[i] - natPdb[j].ca.z));
			
			if (dca_4 <= 4.0) {
                                        d_GDT4 = d_GDT4 + 1.0;
			}
			
			if (dca_2 <= 2.0) {
                                        d_GDT2 = d_GDT2 + 1.0;
                        }
			
			if (dca_1 <= 1.0) {
                                        d_GDT1 = d_GDT1 + 1.0;
                        }
			
			if (dca_05 <= 0.5) {
                                        d_GDT05 = d_GDT05 + 1.0;
                        }

			//GDC-SC calculation using SC atoms
			rsc = sqrt((sc.xt[i] - natPdb[j].sc.x)*(sc.xt[i] - natPdb[j].sc.x) + (sc.yt[i] - natPdb[j].sc.y)*(sc.yt[i] - natPdb[j].sc.y) + (sc.zt[i] - natPdb[j].sc.z)*(sc.zt[i] - natPdb[j].sc.z));
			if (rsc <= 5.0) {
                                r_GDT5 = r_GDT5 + 1.0;
                                if (rsc <= 4.5) {
                                        r_GDT45 = r_GDT45 + 1.0;
                                        if (rsc <= 4.0) {
                                                r_GDT4 = r_GDT4 + 1.0;
                                                if (rsc <= 3.5) {
                                                        r_GDT35 = r_GDT35 + 1.0;
							if (rsc <= 3.0) {
                                                        	r_GDT3 = r_GDT3 + 1.0;
								if (rsc <= 2.5) {
                                                        		r_GDT25 = r_GDT25 + 1.0;
									if (rsc <= 2.0) {
                                                        			r_GDT2 = r_GDT2 + 1.0;
										if (rsc <= 1.5) {
                                                        				r_GDT15 = r_GDT15 + 1.0;
                       									if (rsc <= 1.0) {
                                                 			       			r_GDT1 = r_GDT1 + 1.0;
                                                						if (rsc <= 0.5) {
                                                        						r_GDT05 = r_GDT05 + 1.0;
												}
											}
										}
									}
								}
							}
                                        	}
                                	}
                        	}
			}
			
			point3d p1,p2,p3,p4;
                	p1.x = ca.xt[i];
                	p1.y = ca.yt[i];
                	p1.z = ca.zt[i];
                	p2.x = sc.xt[i];
                	p2.y = sc.yt[i];
                	p2.z = sc.zt[i];
                	p3.x = natPdb[j].ca.x;
                	p3.y = natPdb[j].ca.y;
                	p3.z = natPdb[j].ca.z;
                	p4.x = natPdb[j].sc.x;
                	p4.y = natPdb[j].sc.y;
                	p4.z = natPdb[j].sc.z;
		
			//get distance between CA and SC atoms in model	
			point3d dij1 = getDifference(p1, p2);
			//get distance between CA and SC atoms in native
			point3d dij2 = getDifference(p3, p4);
			//get distance between model's SC atom and native's SC atom
			point3d drij = getDifference(p2, p4);
			if (round(drij.x) == 0 && round(drij.y) == 0 && round(drij.z) == 0) {
                                drij.x = round(drij.x);
                                drij.y = round(drij.y);
                                drij.z = round(drij.z);
                        }

			//get unit vectors i.e. the directions of CA and SC atoms
			point3d uij1 = getUnit(dij1);
			point3d uij2 = getUnit(dij2);
			point3d urij2 = getUnit(drij);
			double omegaij1 = getDotProduct(uij1, urij2);
                	double omegaij2 = getDotProduct(uij2, urij2);
			double omegaij12 = getDotProduct(uij1, uij2);
			
			//get the planar angles
			double thetaij1_rad = acos(omegaij1);
			double thetaij2_rad = acos(omegaij2);
			double thetaij1 = thetaij1_rad * (180/3.1415);
			double thetaij2 = thetaij2_rad * (180/3.1415);

			//get the dihedral angle
			double cosphiij = (omegaij12 - (omegaij1 * omegaij2))/(sin(thetaij1_rad) * sin(thetaij2_rad));
			double phiij_rad = acos(float(cosphiij));
			double phiij = phiij_rad * (180/3.1415);
			phiij = phiij + 180;
			
			//compute planar angle threshold counts
			thetaij1 = round(thetaij1 * 100.0)/100.0;
			if (thetaij1 <= 120) {
				n_thetaij1_120 = n_thetaij1_120 + 1.0;
				if (thetaij1 <= 90) {
					n_thetaij1_90 = n_thetaij1_90 + 1.0;
					if (thetaij1 <= 60) {
						n_thetaij1_60 = n_thetaij1_60 + 1.0;
						if (thetaij1 <= 30) {
							n_thetaij1_30 = n_thetaij1_30 + 1.0;
						}
					}
				}					
			}
			
			thetaij2 = round(thetaij2 * 100.0)/100.0;	
			if (thetaij2 <= 120) {
                                n_thetaij2_120 = n_thetaij2_120 + 1.0;
                                if (thetaij2 <= 90) {
                                        n_thetaij2_90 = n_thetaij2_90 + 1.0;
                                        if (thetaij2 <= 60) {
                                                n_thetaij2_60 = n_thetaij2_60 + 1.0;
                                                if (thetaij2 <= 30) {
                                                        n_thetaij2_30 = n_thetaij2_30 + 1.0;
                                                }
                                        }
                                }
                        }
			
			//compute the dihedral angle threshold counts
			if (phiij <= 300) {
                                n_phiij_300 = n_phiij_300 + 1.0;
                                if (phiij <= 270) {
                                        n_phiij_270 = n_phiij_270 + 1.0;
                                        if (phiij <= 240) {
                                                n_phiij_240 = n_phiij_240 + 1.0;
                                                if (phiij <= 210) {
                                                        n_phiij_210 = n_phiij_210 + 1.0;
							if (phiij <= 180) {
                                                        	n_phiij_180 = n_phiij_180 + 1.0;
								if (phiij <= 150) {
                                                        		n_phiij_150 = n_phiij_150 + 1.0;
									if (phiij <= 120) {
                                                                        	n_phiij_120 = n_phiij_120 + 1.0;
										if (phiij <= 90) {
                                                                        		n_phiij_90 = n_phiij_90 + 1.0;
											if (phiij <= 60) {
                                                                        			n_phiij_60 = n_phiij_60 + 1.0;
												if (phiij <= 30) {
                                                                        				n_phiij_30 = n_phiij_30 + 1.0;
												}
											}
										}
									}
								}
							}
                                                }
                                        }
                                }
                        }
	}		
	
	//calculate the CA component of SPECS
	d_GDT_HA = (d_GDT05 + d_GDT1 + d_GDT2 + d_GDT4)/float(4 * natPdb.size());

	//calculate the SC component of SPECS
	r_GDT05 = r_GDT05/float(natPdb.size());
	r_GDT1 = r_GDT1/float(natPdb.size());
	r_GDT15 = r_GDT15/float(natPdb.size());
	r_GDT2 = r_GDT2/float(natPdb.size());
	r_GDT25 = r_GDT25/float(natPdb.size());
	r_GDT3 = r_GDT3/float(natPdb.size());
	r_GDT35 = r_GDT35/float(natPdb.size());
	r_GDT4 = r_GDT4/float(natPdb.size());
	r_GDT45 = r_GDT45/float(natPdb.size());
	r_GDT5 = r_GDT5/float(natPdb.size());
	
	r_GDC_SC = (1 * 2 * ((10 * r_GDT05) + (9 * r_GDT1) + (8 * r_GDT15) + (7 * r_GDT2) + (6 * r_GDT25) + (5 * r_GDT3) + (4 * r_GDT35) + (3 * r_GDT4) + (2 * r_GDT45) + (1 * r_GDT5)))/(float(11 * 10));

	//calculate theta_1
	n_thetaij1_120 = n_thetaij1_120/float(natPdb.size());
	n_thetaij1_90 = n_thetaij1_90/float(natPdb.size());
	n_thetaij1_60 = n_thetaij1_60/float(natPdb.size());
	n_thetaij1_30 = n_thetaij1_30/float(natPdb.size());
	
	theta_1 = 0.0;	
	if (n_thetaij1_120 == 1 && n_thetaij1_90 == 1 && n_thetaij1_60 == 0 && n_thetaij1_30 == 0) {
		theta_1 = (1 * 2 * ((2 * n_thetaij1_90) + (1 * n_thetaij1_120)))/(float(3 * 2));
	}
	else {
		theta_1 = (1 * 2 * ((4 * n_thetaij1_30) + (3 * n_thetaij1_60) + (2 * n_thetaij1_90) + (1 * n_thetaij1_120)))/(float(5 * 4));
	}

	//calculate theta_2
	n_thetaij2_120 = n_thetaij2_120/float(natPdb.size());
        n_thetaij2_90 = n_thetaij2_90/float(natPdb.size());
        n_thetaij2_60 = n_thetaij2_60/float(natPdb.size());
        n_thetaij2_30 = n_thetaij2_30/float(natPdb.size());
	
	theta_2 = 0.0;
	if (n_thetaij2_120 == 1 && n_thetaij2_90 == 1 && n_thetaij2_60 == 0 && n_thetaij2_30 == 0) {
		theta_2 = (1 * 2 * ((2 * n_thetaij2_90) + (1 * n_thetaij2_120)))/(float(3 * 2));
	}
	else {
        	theta_2 = (1 * 2 * ((4 * n_thetaij2_30) + (3 * n_thetaij2_60) + (2 * n_thetaij2_90) + (1 * n_thetaij2_120)))/(float(5 * 4));
	}

	//calculate phi
        n_phiij_300 = n_phiij_300/float(natPdb.size());
	n_phiij_270 = n_phiij_270/float(natPdb.size());
	n_phiij_240 = n_phiij_240/float(natPdb.size());
        n_phiij_210 = n_phiij_210/float(natPdb.size());
	n_phiij_180 = n_phiij_180/float(natPdb.size());
        n_phiij_150 = n_phiij_150/float(natPdb.size());
	n_phiij_120 = n_phiij_120/float(natPdb.size());
        n_phiij_90 = n_phiij_90/float(natPdb.size());
	n_phiij_60 = n_phiij_60/float(natPdb.size());
        n_phiij_30 = n_phiij_30/float(natPdb.size());
	
	phi = 0.0;
	if (n_phiij_300 == 1 && n_phiij_270 == 1 && n_phiij_240 == 1 && n_phiij_210 == 1 && n_phiij_180 == 1 && n_phiij_150 == 0 && n_phiij_120 == 0 && n_phiij_90 == 0 && n_phiij_60 == 0 && n_phiij_30 == 0) {
		 phi = (1 * 2 * ((5 * n_phiij_180) + (4 * n_phiij_210) + (3 * n_phiij_240) + (2 * n_phiij_270) + (1 * n_phiij_300)))/(float(6 * 5));
	}
	else{
		phi = (1 * 2 * ((10 * n_phiij_30) + (9 * n_phiij_60) + (8 * n_phiij_90) + (7 * n_phiij_120) + (6 * n_phiij_150) + (5 * n_phiij_180) + (4 * n_phiij_210) + (3 * n_phiij_240) + (2 * n_phiij_270) + (1 * n_phiij_300)))/(float(11 * 10));
	}

	SPECS = (4 * d_GDT_HA + r_GDC_SC + theta_1 + theta_2 + phi)/8.0;
        
}

/*************************************************************************
 *  Name        : score_fun
 *  Purpose     : collect those residue with dis<d and calculate scores
 *  Arguments   : vector<pdbInfo> &xt, vector<pdbInfo> &natPdb
 *  Return Type : void
 *************************************************************************/
void score_fun (struc s,vector<pdbInfo> &natPdb) {
	double d_tmp = d;
	long double score_maxsub_sum = 0.0; //MaxSub Score
	long double score_sum = 0.0;	//TMScore
	long double score_sum10 = 0.0; //TMScore10
	long double dis = 0.0;	
	
	for (;;) {	
		n_cut=0;                   //number of residue-pairs dis<d, for iteration
      		n_GDT05=0.0;                 //for GDT-score, # of dis<0.5
      		n_GDT1=0.0;                  //for GDT-score, # of dis<1
      		n_GDT2=0.0;                  //for GDT-score, # of dis<2
     		n_GDT4=0.0;                  //for GDT-score, # of dis<4
      		n_GDT8=0.0;                  //for GDT-score, # of dis<8
      		score_maxsub_sum = 0.0;
		score_sum = 0.0;
		score_sum10 = 0.0;

		for (int k = 1; k <= align.n_ali; k++) {	
			int i = align.iA[k]; // ![1,nseqA] reorder number of structure A
			int j = align.iB[k]; // ![1,nseqB]
			dis = sqrt((s.xt[i] - natPdb[j].ca.x)*(s.xt[i] - natPdb[j].ca.x) + (s.yt[i] - natPdb[j].ca.y)*(s.yt[i] - natPdb[j].ca.y) + (s.zt[i] - natPdb[j].ca.z)*(s.zt[i] - natPdb[j].ca.z)); 
			if (dis < d_tmp) {	
				n_cut = n_cut + 1;
				i_ali[n_cut] = k;  // ![1,n_ali], mark the residue-pairs in dis<d
			}
			//for calculating GDT-score thresholds
			if (dis <= 8.0) {
            			n_GDT8 = n_GDT8 + 1.0;
            			if (dis <= 4.0) {
               				n_GDT4 = n_GDT4 + 1.0;
               				if (dis <= 2.0) {
                  				n_GDT2 = n_GDT2 + 1.0;
                  				if (dis <= 1.0) {
                     					n_GDT1 = n_GDT1 + 1.0;
                     					if(dis <= 0.5) {
                        					n_GDT05 = n_GDT05 + 1.0;
                     					}
						}
					}
				}
			}
			//for calculating MaxSub-score
			if (dis < 3.5) {
				score_maxsub_sum = score_maxsub_sum + 1.0/(1.0+((dis/3.5)*(dis/3.5)));
			}	
			//for calculating TMScore
			score_sum = score_sum + 1.0/(1.0 + (dis/d0)*(dis/d0));
			
			//for TMScore10
			if (dis < 10.0) {
				score_sum10 = score_sum10 + 1.0/(1.0+(dis/d0)*(dis/d0));
			}
		}
	
		if ((n_cut >= 3) || (align.n_ali <= 3)) {
			break;
		}

		d_tmp = d_tmp + 0.5;
		continue;
		break;
	}				
	//TMScore
	score = score_sum/float(natPdb.size());
	//MaxSub Score
	score_maxsub = score_maxsub_sum/float(natPdb.size());
	//TMScore10
	score10 = score_sum10/float(natPdb.size());
}

/*************************************************************************
 * Name        : u3b
 * Purpose     : calculate sum of (r_d-r_m)^2
 * 		 w is weight for atom pair c m
 * 		 x(i,m) are coordinates of atom c m in set x
 * 		 y(i,m) are coordinates of atom c m in set y
 * 		 n is number of atom pairs
 * 		 mode-0: calculate rms only
 * 		 mode-1: calculate rms, u, t
 * 		 rms-sum of w*(ux*t-y)**2 over all atom pairs
 * 		 u-u(i,j) is rotation matrix for best superposition
 * 		 t-t(i) is translation matrix for best superposition
 * 		 ier-0:a unique optimal superposition has been determined
 * 		     1:superposition is not unique but optimal
 * 		     2:no result obtained because of negative weights w or all weights equal to zero
 * Arguments   : double w[3000], double x[4][3000], double y[4][3000], int n, int mode, double rms, double u[4][4], double t[4], int ier
 * Return Type : void
 *************************************************************************/
long double u3b(long double w[3000], long double x[4][3000], long double y[4][3000], int n, int mode, long double rms, long double u[4][4], long double t[4], int ier) {
	int m1,j;
	long double sqrt3 = 1.73205080756888;
        long double tol = 0.01; 
	long double zero = 0.0;
	int ip[] = {-100, 1, 2, 4, 2, 3, 5, 4, 5, 6};
	int ip2312[] = {-100, 2, 3, 1, 2};
	long double r[3][3], xc[3], yc[3], wc;
	long double a[3][3], b[3][3], rr[6], ss[6];
	long double e0, d, det, h, g;
	long double e[3];
	long double cth, sth, sqrth, p, sigma;
        long double c1x, c1y, c1z, c2x, c2y, c2z;
	long double s1x, s1y, s1z, s2x, s2y, s2z;
	long double sxx, sxy, sxz, syx, syy, syz, szx, szy, szz;		

	//initialization	
	wc = zero;
	rms = zero;
	e0 = zero;
	s1x = zero;
      	s1y = zero;
      	s1z = zero;
      	s2x = zero;
      	s2y = zero;
      	s2z = zero;
      	sxx = zero;
      	sxy = zero;
      	sxz = zero;
      	syx = zero;
      	syy = zero;
      	syz = zero;
      	szx = zero;
      	szy = zero;
      	szz = zero;

	for (int i = 1; i <= 3; i++) {
		xc[i] = zero;
		yc[i] = zero;
		t[i] = zero;
		for (int j = 1; j <= 3; j++) {
			u[i][j] = zero;
			r[i][j] = zero;
			a[i][j] = zero;
			if (i == j) {
				u[i][j] = 1.0;
				a[i][j] = 1.0;
			}
		}
	}
	
	ier = -1;
	if (n < 1) {
		return 0.0;
	}

	ier = -2;
	for (int m = 1; m <= n; m++) {
		c1x=x[1][m];
         	c1y=x[2][m];
         	c1z=x[3][m];

         	c2x=y[1][m];
         	c2y=y[2][m];
         	c2z=y[3][m];

         	s1x = s1x + c1x;
         	s1y = s1y + c1y;
         	s1z = s1z + c1z;

         	s2x = s2x + c2x;
         	s2y = s2y + c2y;
         	s2z = s2z + c2z;

         	sxx = sxx + c1x*c2x;
         	sxy = sxy + c1x*c2y;
         	sxz = sxz + c1x*c2z;

         	syx = syx + c1y*c2x;
         	syy = syy + c1y*c2y;
         	syz = syz + c1y*c2z;

         	szx = szx + c1z*c2x;
         	szy = szy + c1z*c2y;
         	szz = szz + c1z*c2z;
	}
	
	xc[1] = s1x/n;
      	xc[2] = s1y/n;
      	xc[3] = s1z/n;

      	yc[1] = s2x/n;
      	yc[2] = s2y/n;
      	yc[3] = s2z/n;
	
	if (mode == 2 || mode == 0) {
		for (int m = 1; m <= n; m++) {
			for (int i = 1; i <= 3; i++) {
				e0 = e0+ (x[i][m]-xc[i])*(x[i][m]-xc[i]) + (y[i][m]-yc[i])*(y[i][m]-yc[i]);
			}
		}
	}
	
	r[1][1] = sxx-s1x*s2x/n;
	r[2][1] = sxy-s1x*s2y/n;
 	r[3][1] = sxz-s1x*s2z/n;
	r[1][2] = syx-s1y*s2x/n;
        r[2][2] = syy-s1y*s2y/n;
        r[3][2] = syz-s1y*s2z/n;
	r[1][3] = szx-s1z*s2x/n;
        r[2][3] = szy-s1z*s2y/n;
        r[3][3] = szz-s1z*s2z/n;
	
	det = r[1][1] * (r[2][2] * r[3][3] - r[2][3] * r[3][2]) - r[1][2] * (r[2][1] * r[3][3] - r[2][3] * r[3][1]) + r[1][3] * (r[2][1] * r[3][2] - r[2][2] * r[3][1]);
	sigma = det;
	
	//compute tras(r)*r
	int m = 0;
	for (j = 1; j <= 3; j++) {
		for (int i = 1; i <= j; i++) {
			m = m + 1;
			rr[m] = (r[1][i] * r[1][j] + r[2][i] * r[2][j] + r[3][i] * r[3][j]);
		}
	}
	double spur = (rr[1] + rr[3] + rr[6]) / 3.0;
	double cof = (((((rr[3]*rr[6] - rr[5]*rr[5]) + rr[1]*rr[6]) - rr[4]*rr[4]) + rr[1]*rr[3]) - rr[2]*rr[2]) / 3.0;
	det = det * det;
	for (int i = 1; i <= 3; i++) {
		e[i] = spur;
	}

	if (spur > zero) {
		d = spur * spur;
		h = d - cof;
		g = (spur * cof - det) / 2.0 - spur * h;
		if (h > zero) {
			sqrth = sqrt(h);
			d = h * h * h - g * g;
			if (d < zero) {
				d = zero;
			}
			
			d = atan2(sqrt(d),-g)/3.0;
			cth = sqrth * cos(d);
			sth = sqrth * sqrt3 * sin(d);
			e[1] = (spur + cth + cth);
			e[2] = (spur - cth + sth);
			e[3] = (spur - cth - sth);
			if (mode != 0) { //compute a
				for (int l = 1; l <= 3; l = l + 2) {
					d = e[l];
					ss[1] = ((d - rr[3]) * (d - rr[6]) - rr[5] * rr[5]);
					ss[2] = ((d - rr[6]) * rr[2] + rr[4] * rr[5]);
					ss[3] = ((d - rr[1]) * (d - rr[6]) - rr[4] * rr[4]);
					ss[4] = ((d - rr[3]) * rr[4] + rr[2] * rr[5]);
					ss[5] = ((d - rr[1]) * rr[5] + rr[2] * rr[4]);
					ss[6] = ((d - rr[1]) * (d - rr[3]) - rr[2] * rr[2]);
						
					if (fabs(ss[1]) >= fabs(ss[3])) {
						j = 1;
						if (fabs(ss[1]) < fabs(ss[6])) {
							j = 3;
						}
					}
					else if (fabs(ss[3]) >= fabs(ss[6])) {
						j = 2;
					}
					else {
						j = 3;
					}
					
					d = 0.0;
					j = 3 * (j - 1);
					for (int i = 1; i <= 3; i++) {
						int k = ip[i + j];
						a[i][l] = ss[k];
						d = d + ss[k] * ss[k];
					}
					
					if (d > 0.0) {
						d = 1.0/sqrt(d);
					} 
					else {
						d = 0.0;
					}
					
					for (int i = 1; i <= 3; i++) {
						a[i][l] = a[i][l] * d;
					}
				}//for l
				
				d = ((a[1][1] * a[1][3]) + (a[2][1] * a[2][3]) + (a[3][1] * a[3][3]));
				if ((e[1] - e[2]) > (e[2] - e[3])) {
					m1 = 3;
					m = 1;
				}
				else {
					m1 = 1;
					m = 3;
				}
				
				p = zero;
				for (int i = 1; i <= 3; i++) {
					a[i][m1] = a[i][m1] - d * a[i][m];
					p = p + a[i][m1] * a[i][m1];
				}
				
				if (p <= tol) {
					p = 1.0;
					for (int i = 1; i <= 3; i++) {
						if (p >= fabs(a[i][m])) {
							p = fabs(a[i][m]);
							j = i;
						}
					}
					
					int k = ip2312[j];
					int l = ip2312[(j + 1)];
					p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
					if (p > tol) {
						a[j][m1] = 0.0;
						a[k][m1] = (-a[l][m]/p);
						a[k][m] = a[k][m]/p;
					}
				}
				else {
					p = 1.0/sqrt(p);
					for (int i = 1; i <= 3; i++) {
						a[i][m1] = a[i][m1] * p;
					}
				}
				a[1][2] = a[2][3] * a[3][1] - a[2][1] * a[3][3];
				a[2][2] = a[3][3] * a[1][1] - a[3][1] * a[1][3];
				a[3][2] = a[1][3] * a[2][1] - a[1][1] * a[2][3];
			}//if(mode != 0)
			
			
		}//h>0
		
		if (mode != 0) {	
			//compute b anyway
			for (int l = 1; l <= 2; l++) {
				d = 0.0;
				for (int i = 1; i <= 3; i++) {
					b[i][l] = (r[i][1] * a[1][l] + r[i][2] * a[2][l] + r[i][3] * a[3][l]);
					d = d + b[i][l] * b[i][l];
				}
				if (d > 0.0) {
					d = 1.0/sqrt(d);
				}
				else {
					d = 0.0;
				}	
			
				for (int i = 1; i <= 3; i++) {
					b[i][l] = b[i][l] * d;
				}
			}
		
			d = b[1][1] * b[1][2] + b[2][1] * b[2][2] + b[3][1] * b[3][2];
			p = 0.0;
		
			for (int i = 1; i <= 3; i++) {
				b[i][2] = b[i][2] - d * b[i][1];
				p = p + b[i][2] * b[i][2];
			}
			
			if (p <= tol) {
				p = 1.0;
				for (int i = 1; i <= 3; i++) {
					if (p >= fabs(b[i][1])) {
						p = fabs(b[i][1]);
						j = i;
					}
				}

				int k = ip2312[j];
				int l = ip2312[j + 1];
				p = sqrt(b[k][1] * b[k][1] + b[l][1] * b[l][1]);
				if (p > tol) {
					b[j][2] = 0.0;
					b[k][2] = (-b[l][1] / p);
					b[l][2] = (b[k][1] / p);
				}	
			}				
			else {
				p = 1.0/sqrt(p);
				for (int i = 1; i <= 3; i++) {
					b[i][2] = b[i][2] * p;
				}
			}
		}
		
		b[1][3] = b[2][1] * b[3][2] - b[2][2] * b[3][1];
		b[2][3] = b[3][1] * b[1][2] - b[3][2] * b[1][1];
		b[3][3] = b[1][1] * b[2][2] - b[1][2] * b[2][1];
		//compute u
		for (int i = 1; i <= 3; i++) {
			for (int j = 1; j <= 3; j++) {
				u[i][j] = b[i][1] * a[j][1] + b[i][2] * a[j][2] + b[i][3] * a[j][3];
			}
		}
		
		//compute t
		for (int i = 1; i <= 3; i++) {
                        t[i] = (yc[i] - u[i][1] * xc[1] - u[i][2] * xc[2] - u[i][3] * xc[3]);
                }
	}
	else {
		for (int i = 1; i <= 3; i++) {
			t[i] = (yc[i] - u[i][1] * xc[1] - u[i][2] * xc[2] - u[i][3] * xc[3]);
		}
	}//else spur > 0
	
	//compute rms
	for (int i = 1; i <= 3; i++) {
		if (e[i] < 0) {
			e[i] = 0;
		}
		e[i] = sqrt(e[i]);
	}
	
	ier = 0;
	if (e[2] <= (e[1] * 1.0 - 0.5)) 
		ier = -1;
	
	d = e[3];
	if (sigma < 0.0) {
		d = -d;
		if ((e[2] - e[3]) <= (e[1] * 1.0 - 0.5)) {
			ier = -1;
		}
	}
	
	d = (d + e[2]) + e[1];
	if (mode == 2 || mode == 0) {
		rms = (e0 - d) - d;
		if (rms < 0.0)
			rms = 0.0;
		return rms;
	}
}

/*************************************************************************
 * Name        : getDotProduct
 * Purpose     : gets the dot product (i.e. scalar product) for two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : double
 *************************************************************************/
double getDotProduct(point3d & p1, point3d &p2) {
    return (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
}

/*************************************************************************
 * Name        : getNorm
 * Purpose     : gets the norm (i.e. length) of a vector
 * Arguments   : point3d & p
 * Return Type : double
 *************************************************************************/
double getNorm(point3d & p) {
    return sqrt( p.x * p.x + p.y * p.y + p.z * p.z);
}

/*************************************************************************
 * Name        : getDifference
 * Purpose     : gets the difference vector between two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : point3d
 *************************************************************************/
point3d getDifference(point3d & p1, point3d &p2) {
    point3d p;
    p.x = p2.x - p1.x;
    p.y = p2.y - p1.y;
    p.z = p2.z - p1.z;
    return p;
}

/*************************************************************************
 * Name        : getUnit
 * Purpose     : gets the unit vector for a vector
 * Arguments   : point3d & p
 * Return Type : point3d
 *************************************************************************/
point3d getUnit(point3d & p) {
    point3d u;
    double norm = getNorm(p);
    if (norm == 0.0) {
        u.x = 0.0;
        u.y = 0.0;
        u.z = 0.0;
    }
    else {
        u.x = p.x / norm;
        u.y = p.y / norm;
        u.z = p.z / norm;
    }
    return u;
}

/*************************************************************************
 * Name        : getCentroid
 * Purpose     : gets centroid of a cloud of points
 * Arguments   : vector<point3d> &pointCloud
 * Return Type : point3d
 *************************************************************************/
point3d getCentroid(vector<point3d> &pointCloud) {
    point3d centroid;
    centroid.x = 0.0;
    centroid.y = 0.0;
    centroid.z = 0.0;
    for (int i = 0; i < pointCloud.size(); i++) {
        centroid.x += pointCloud[i].x / pointCloud.size();
        centroid.y += pointCloud[i].y / pointCloud.size();
        centroid.z += pointCloud[i].z / pointCloud.size();
    }
    return centroid;
}
