
/*
 *
 *  Created by Hilal Kosucu on 15/01/09.
 *
 *  Extended by Yaron Orenstein on 15/01/16.
 *  
 *  Extended by Shu Yang on 15/04/17.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "lib/lbfgsb.h"
#include "rnacontextLibrary.cpp"
#include <iomanip>
#include <assert.h>



using namespace std;


//
int iter = 0;
static ap::real_1d_array AVAL;

//prototypes
double calcScore(int seqID);
void findBestScores(string filename, int motifWidth, int numBest);
double calcScoreTest(int seqID);
double calcScoreTestWindow(int seqID, int window);
double calcScoreTestSeq(int seqID);
int readData(string inputFileName, vector<string> & seqs, vector<double> & ratios);
int readModel(string inputFileName, double *baseParams, double **annotParams, double & biasS, double & biasA, double & alpha, double & beta, int mw, int size);

// Helper function: converts bases e.g. A,C,G,U to numbers 0,1,2,3
int basetoInt (char base)
{
	for(int i = 0 ; i < AlphabetSize ; i++)
   	{ 
   		if(base == baseAlphabet[i])
			return  i;
	}
	return -1;   // Error 
	
}

char inttoBase(int i)
{
	return baseAlphabet[i];
}

int getInt(const char *str, int k, int motifWidth) {
        int rc=0;
        for (int n = 0; n < motifWidth; n++)
                rc = rc + pow(4, n) * basetoInt(str[k+n]);
        return rc;
}

string getString(int i, int motifWidth) {
	string rc="";
	for (int j = 0; j < motifWidth; j++) {
		rc = rc + inttoBase(i % 4);
		i = i / 4;
	}
	return rc;
}

// convert PWM value (from expwm[][]; 1st dimension is ACGU, 2nd dimension is position) to k-mer value (i.e. entry in baseParams[])
double getPWMdouble(double **pwm, int i, int motifWidth){
	double sum = 0.0;

	for (int j = 0; j < motifWidth; j++) {
		sum += pwm[i % 4][j];
		i = i / 4;
	}
	return sum;
}

// Helper function: converts annotation alphabet to numeric values 
int annottoInt (char annotAlph)
{
	for(int i = 0 ; i < AnnotAlphabetSize ; i++)
   	{ 
   		if(annotAlph == annotAlphabet[i])
			return  i;
	}
	return -1;   // Error 
}

// Helper Function: sum the corresponding base parameters of a kmer "motif"
double sumBaseParamsFunc(string motif)
{
	return baseParams[getInt(motif.c_str(), 0, motifWidth)];
}

// Helper Function: sum the correponding kmer starting at startIndex within sequence specified by seqID
double sumBaseParamsFunc(int seqID, int startIndex)
{
	double sum = 0;
	return baseParams[getInt(seqsTrain[seqID].c_str(), startIndex, motifWidth)];
}

// Helper Function: has the same functionality as sumBaseParamsFunc but uses test sequences
double sumBaseParamsTestFunc(int seqID, int startIndex)
{
	double sum = 0;
	return baseParams[getInt(seqsTest[seqID].c_str(), startIndex, motifWidth)];
}

// Helper Function: sums the annotation or structure parameters for the kmer starting at "startIndex" in sequence "seqID" for annotation "letter"
double sumProfilesFunc(int seqID, int letter, int startIndex)
{
	double sum = 0;
	for(int i =0; i <motifWidth; i++)
	{    
		 sum += inputAnnotsTrain[seqID][letter][startIndex +i];
	}
	return sum;	
}

// Helper Function: has the same functionality as sumProfilesFunc, but for test sequences
double sumProfilesFuncTest(int seqID, int letter, int startIndex)
{
	double sum = 0;
	for(int i =0; i <motifWidth; i++)
	{    
		 sum += inputAnnotsTest[seqID][letter][startIndex +i];
	}
	return sum;	
}


// Takes the updated parameters from paraVec and restore those values to the structures defined for parameters  
void vecToPara(ap::real_1d_array& paraVec) {
	int count = 1;

	// map back the base parameters
	int pownum = pow(AlphabetSize, motifWidth);
	// comment out below
	//for (int i=0; i<pownum; i++){
	//  baseParams[i] = paraVec(count);
	//  count++;
 //  	} 

   // map back the annot parameters
    for (int i = 0; i < pownum; i++)
	    for (int j=0; j<AnnotAlphabetSize; j++){
		  annotParams[i][j] = paraVec(count);  
		  count++;
    }
	//map back bias
	biasS = paraVec(count); 		
	count++;
	
	biasA = paraVec(count);	
	count++;  
	  
	parA = exp(paraVec(count)); 	
	count++;
	
	parB = paraVec(count); 	
}

// Reads the annotation profile from supplied file with filename "filename"
void readAnnots(char* filename, Annots & annots, int seqNum, int ignoreAnnot = 0, bool meanvar = false)
{
	string line;
	string structure;
	string seq;
	double prob;
	int seqCount = 0;
	int pos = 0;
	int seqLengths[seqNum];
    int sumLengths = 0;
	double data;
    ifstream input;
    input.open(filename);
    if(input.fail())
		cout<<"cannot open the annotation file"<<endl;
	while(getline(input, line))
	{	
		for(int letter = 0 ; letter < AnnotAlphabetSize; letter++)
		{			
			getline(input,line);
			istringstream ins;
			ins.str(line);
			while(!ins.eof())
			{
				ins>>data;
				if (ignoreAnnot == 0) {
					annots[seqCount][letter][pos] = data;
				}
				else {annots[seqCount][letter][pos] = 1 / AnnotAlphabetSize;}
				pos++;
			}
			seqLengths[seqCount] = pos - 1;
			pos = 0;			         	 
		}  	    
		seqCount++;	        
	}
	if(meanvar)
	{
		double sum= 0;
		double diffSum = 0; 
		double diff = 0;
		
		for(int letter = 0; letter < AnnotAlphabetSize; letter++)
		{
			sumLengths = 0;
			for(int seqID = 0; seqID < seqNum ; seqID++)
			{
				for(int pos = 0 ; pos < seqLengths[seqID]; pos++)
				{ 
			    	sum += annots[seqID][letter][pos];

				}
				sumLengths += seqLengths[seqID];
			}			
			means[letter] = sum / sumLengths;
			for(int seqID = 0; seqID < seqNum ; seqID++)
			{
				for(int pos = 0 ; pos < seqLengths[seqID]; pos++)
				{
					diff = means[letter] - annots[seqID][letter][pos];
					diffSum += diff * diff;
				}			
			}
			variances[letter] = diffSum / sumLengths ;
			diffSum = 0;
			sum = 0;
		}
	}
}


// Get the updated parameters from optimization module and recalculate the scores of the sequences
double originalFunc(ap::real_1d_array& paraVec, bool last = false) 
{
	double result = 0;
	double regVal = 0;
	vecToPara(paraVec); // get the updated parameters	// note vecToPara() is modified since paraVec doesnot contain base parameters anymore
	for(int i = 0; i < trainSeqNum; i++)
	{
		scoresTrain[i] = calcScore(i);	// note calcScore() is implemented differently from the paper (both RCK and RNAcontext papers), see comments inside that function
	}
	for (int i = 0; i < trainSeqNum; i++) 
	{
		result += square(ratiosTrain[i] - parA * scoresTrain[i] - parB) ;	// this is the loss term (first square term) in formula 5 of RCK paper
	}

	// regularization
	int pownum= pow(AlphabetSize, motifWidth);
	for (int l=0; l<pownum; l++){
			regVal += 1 / (1+exp(-biasS-baseParams[l]));	// although baseParams[] is a constant now, I still keep the regularization term the same as RCK (which is different from RNAcontext)
	}
	for (int l = 0; l < pownum; l++)
		for (int k = 0; k < AnnotAlphabetSize; k++) {	// note the RCK implementation leaves this curly bracket empty which may indicate they have tried to regularize annotParams[][] as well but decided not to do it.
		}
     
	result += regVal * SMALL2;   // SMALL2 corresponds to alpha
	return result;
}

// Test the parameters learned by the optimization procedure by scpring test sequences and output these scores
void testBestParams(string outFileName, string modelFileName, int numParams)
{
	ofstream out;
	out.open(outFileName.c_str());
	ofstream modelOut;
	modelOut.open(modelFileName.c_str());
	double result = 0;
	for(int i = 0; i < testSeqNum ; i++)
	{ 	   
   	 	scoresTest[i] = calcScoreTest(i);  
     	out<<ratiosTest[i]<<"\t"<<scoresTest[i]<<endl;  	
	}  

	for (int i = 0; i < testSeqNum; i++) 
	{
		result += square(ratiosTest[i] - parA * scoresTest[i] - parB) ;
	}
	modelOut<<"Error on the test set: "<<result<<endl;
}
// Test the parameters learned by the optimization procedure by scpring test sequences and output these scores
void testBestParamsSeq(string outFileName, string modelFileName, int numParams)
{
	double score = 0;
	double result = 0;	
	ofstream out;
	out.open(outFileName.c_str());
	ofstream modelOut;
	modelOut.open(modelFileName.c_str());

	
	for(int i = 0; i < testSeqNum ; i++)
	{ 	   
   	 	score = calcScoreTestSeq(i);  
	        out<<ratiosTest[i]<<"\t"<<score<<endl; 
		result += square(ratiosTest[i] - parA * score - parB) ;
	}  

	modelOut<<"Error on the test set seq only: "<<result<<endl;
}

// Initializes and optimized the parameters	// pass in additional external pwm file name as parameter to the function
void searchMotif(int seedNum, int seedWidth, int minWidth, int maxWidth, char *outFileName, char* annotFileName, char *outDir, int maxiter, int model, char *expwmFileName) 
{
		
	ofstream modelOut;
	ofstream outTrain;
		
	bool once = false;
	string annotFileString = annotFileName;
	string annotAlphStr = annotAlphabet;
	string outFile = outFileName;
	string dir = outDir;
	string modelFileName;
	string testFileName;
	string testFileNameSeq;
	string rankFileName; 
	int numParams;

	int initNum;
	int inisave;
	
	double smallFactor = 0.1;
	double error = 0;
	double bestError = MAX_ERROR;
	
	int numBest = 4;
	
	int alphIndex, annotIndex;
	string paramFilename;
	paramFilename = dir + "/params_" + outFile + ".txt" ;
	ofstream outParam;
	outParam.open(paramFilename.c_str());
	if(outParam.fail())
		cout<<paramFilename<<" can not be opened "<<endl;
	/**************** optimizer settings **********************************/ 
	const double epsg     = 1e-6;
	const double epsf     = 1e-6;
	const double epsx     = 1e-6;
//	const int maxiter  = 200;
	const int m = 5;
	int info; // holds the information about why search is terminated
	int paramCounter;
	int numDim;  //number of dimensions
	/**************** optimizer settings **********************************/ 
	// used to record the results with smaller motifWidth.
	double biasAPrev, biasSPrev, parAPrev, parBPrev;
	//best parameters
	double bestBiasA, bestBiasS, bestParA, bestParB;
	// K-mers
	baseParams = (double *)malloc(pow(AlphabetSize, maxWidth)*sizeof(double));
	annotParams = (double **)malloc(pow(AlphabetSize, maxWidth)*sizeof(double*));
	for (int i = 0; i < pow(AlphabetSize, maxWidth); i++)
		annotParams[i] = (double *)malloc(AnnotAlphabetSize*sizeof(double));

	double *baseParamsPrev; double **annotParamsPrev;
	double *bestBaseParams, *bestBaseParamsScaled; double **bestAnnotParams;
	
	/********ALLOCATE MEMORY FOR PARAMSPREV*****************************************/
	baseParamsPrev = (double *)malloc(pow(AlphabetSize, maxWidth)*sizeof(double));

	// read in external PWMs parameters from pwm file. The PWMs should contain log-transfered value (should be consistent with sumBaseParamsFunc() in RNAcontext (not RCK))
	double **expwm = (double **)malloc(AlphabetSize * sizeof(double*));
	for (int i=0; i< AlphabetSize; i++)
		expwm[i] = (double *)malloc(maxWidth * sizeof(double));	// here, maxWidth is expected to be the same as minWidth
	FILE *expwmfp = fopen(expwmFileName, "r");
	if (expwmfp == NULL) {
		perror(" fail to read pwm file");
		exit(-1);
	}
	ReadPWM(expwmfp, maxWidth, AlphabetSize, expwm);


	/**************** BIG LOOP MOTIFWIDTH CHANGES **********************************/
	for (int mw = minWidth; mw <= maxWidth; mw++)
	{
		motifWidth = mw;
		int pownum = pow(AlphabetSize, motifWidth);
		bestError = MAX_ERROR;

		baseParams = (double *)malloc(pow(AlphabetSize, maxWidth)*sizeof(double));
		bestBaseParams = (double *)malloc(pow(AlphabetSize, maxWidth)*sizeof(double));
		bestBaseParamsScaled = (double *)malloc(pow(AlphabetSize, maxWidth)*sizeof(double));
		bestAnnotParams = (double **)malloc(pow(AlphabetSize, maxWidth) * sizeof(double*));
		annotParamsPrev = (double **)malloc(pow(AlphabetSize, maxWidth) * sizeof(double*));
		for (int i = 0; i < pow(AlphabetSize, maxWidth); i++) {
			bestAnnotParams[i] = (double*)malloc(AnnotAlphabetSize * sizeof(double));
			annotParamsPrev[i] = (double*)malloc(AnnotAlphabetSize * sizeof(double));
		}
		once = false;
		   
		initNum = MAX(3, seedNum);

		/* optimizer settings starts here */

		numDim = pow(AlphabetSize, motifWidth) * (AnnotAlphabetSize) + 2 + 2;// numDim = pow(AlphabetSize, motifWidth) * (1+AnnotAlphabetSize) + 2 + 2;   // +1 for bias +2 for parA and parB
		const int n = numDim;
		ap::real_1d_array x;
		x.setbounds(1,n);
		
		ap::integer_1d_array  constraints;
		ap::real_1d_array lower_bounds;
		ap::real_1d_array upper_bounds;
		
		constraints.setbounds(1,n);
		lower_bounds.setbounds(1,n);
		upper_bounds.setbounds(1,n);
		/* optimizer settings until here  */
		
		/***********Memory Allocation for Test Set***********/
		for(int i = 0 ; i < testSeqNum ; i++)
		{
			kmerNumsTest[i] = seqsTest[i].length()- motifWidth +1; 
		}
		/***********Memory Allocation for Training Set***********/
		for(int i = 0 ; i < trainSeqNum ; i++)
		{	
			kmerNumsTrain[i] = seqsTrain[i].length()- motifWidth +1; 
			affinityValsTrain[i] = (double*) malloc(sizeof(double) * kmerNumsTrain[i]);
			sumBaseParamsTrain[i] = (double*) malloc(sizeof(double) * kmerNumsTrain[i]); 
			AffinitySeqVals[i] = (double*) malloc(sizeof(double) * kmerNumsTrain[i]);
		        AffinityAnnotVals[i] = (double*) malloc(sizeof(double) * kmerNumsTrain[i]);
			sumAnnotParamsTrain[i] = (double**) malloc(sizeof(double*) * AnnotAlphabetSize);
			for(int j =0 ; j < AnnotAlphabetSize ; j++)
			{
				sumAnnotParamsTrain[i][j] = (double*) malloc(sizeof(double) * kmerNumsTrain[i]);     
			}

		}
		
		string mwStr;
		stringstream outs;
		outs << motifWidth;
		mwStr = outs.str();
		/**************** DIFFERENT INITIALIZATIONS **********************************/
		for (int ini = 0; ini < initNum; ini++) 
		{
			paramCounter = 1; //counter for parameter array of the optimizer
			
			if (motifWidth == minWidth) 
			{
				cout<<"MOTIFWIDTH "<<motifWidth<<endl;
				parA = 1;
				parB = RANDOM_PARAM_ENTRY;
				biasS = -2;
				biasA = -2;
				for (int j = 0; j < pownum; j++) {
					baseParams[j] = RANDOM_PARAM_ENTRY; 
				}
				for (int j = 0; j < pownum; j++)
					for (int k = 0; k < AnnotAlphabetSize; k++)
						annotParams[j][k] = RANDOM_PARAM_ENTRY;
			}
			// all the others, based on the previous results
			else 
			{
				inisave= ini % 4;
				switch (inisave) {
					case 0: 
						cout<<inisave << " MOTIFWIDTH "<<motifWidth<<" Initialization "<<ini<<endl;
						biasA = biasAPrev;
						biasS = biasSPrev;
						parA = parAPrev;
						parB = parBPrev;
						for (int j = 0; j < pownum; j++) {
							baseParams[j] = baseParamsPrev[j / 4];
						}
						for (int j = 0; j < pownum; j++)
							for (int k = 0; k < AnnotAlphabetSize; k++)
								annotParams[j][k] = annotParamsPrev[j][k];
						break;

					case 1:
						cout<<inisave << " MOTIFWIDTH "<<motifWidth<<" Initialization "<<ini<<endl;
						
						parA = parAPrev;
						parB = parBPrev;
						biasS = biasSPrev;
						biasA = biasAPrev;
						for (int j = 0; j < pownum; j++) {
							baseParams[j] = baseParamsPrev[j / 4];
						}
						for (int j = 0; j < pownum; j++)
							for (int k = 0; k < AnnotAlphabetSize; k++)
								annotParams[j][k] = annotParamsPrev[j][k];
						break;
						
					case 2:
						cout<<inisave << " MOTIFWIDTH "<<motifWidth<<" Initialization "<<ini<<endl;
						parA = 1;
						parB = RANDOM_PARAM_ENTRY;
						biasS = -2;
						biasA = -2;
						for (int j = 0; j < pownum; j++) {
							baseParams[j] = RANDOM_PARAM_ENTRY * smallFactor;
						}
						for (int j = 0; j < pownum; j++)
							for (int k = 0; k < AnnotAlphabetSize; k++)
								annotParams[j][k] = RANDOM_PARAM_ENTRY*smallFactor;
						break;
					default:
						cout<<inisave << " MOTIFWIDTH "<<motifWidth<<" Initialization "<<ini<<endl;
						
						parA = 1;
						parB = RANDOM_PARAM_ENTRY;
						biasS = -2;
						biasA = RANDOM_PARAM_ENTRY;
						for (int j = 0; j < pownum; j++) {
							baseParams[j] = RANDOM_PARAM_ENTRY * smallFactor;
						}
						for (int j = 0; j < pownum; j++)
							for (int k = 0; k < AnnotAlphabetSize; k++)
								annotParams[j][k] = RANDOM_PARAM_ENTRY * smallFactor;
						break;
				}  //closes switch
			} //closes the else loop which is for non-minMotidWidths
			
			/**************** initialize the parameters of the optimizer **********************************/
			if (model == 1) {	// need to double-check whether/when this condition would be satisfied at all
				modelFileName = dir + "/model_" + outFile + "_" + mwStr + ".txt";
				readModel(modelFileName, baseParams, annotParams, biasS, biasA, parA, parB, pownum, AnnotAlphabetSize);
			}

			// comment out below
			//for (int j=0; j < pownum; j++){ 
			//	x(paramCounter)= baseParams[j]; 	
			//	constraints(paramCounter) = 0; 		
			//	paramCounter++;
			//}
			// assign pwm values to baseParams
			for (int j = 0; j < pownum; j++){
				baseParams[j] = getPWMdouble(expwm, j, motifWidth);
			}

			for (int j = 0; j < pownum; j++)
				for (int k = 0; k < AnnotAlphabetSize; k++) {
					x(paramCounter) = annotParams[j][k];
					constraints(paramCounter) = 0;
					paramCounter++;
				}
			x(paramCounter) =  biasS;
			constraints(paramCounter) = 0; 	
			paramCounter++;
			x(paramCounter) =  biasA;
			constraints(paramCounter) = 0; 	
			paramCounter++;
			x(paramCounter) =  log(parA);
			constraints(paramCounter) = 0;
			paramCounter++;
			x(paramCounter) =  parB;
			constraints(paramCounter) = 0; 	
			paramCounter++;
	
			/****************  LBFGS CALL *****************************************************************/  
			
			lbfgsbminimize(n, m, x, epsg, epsf, epsx, maxiter, constraints, lower_bounds, upper_bounds, info);    
	   		cout<<"The search is terminated. Case "<<info<<endl; // gives info about why the search is terminated
			vecToPara(x);  //copy the resulting parameters into the variables
			error = originalFunc(x);   
			/* update the best likelihood and all the best parameters */
			if (error < bestError) 
			{
				bestError = error;
				bestParA = parA;
				bestParB = parB;
				bestBiasA = biasA;
				bestBiasS = biasS;
				for (int i = 0; i<pownum; i++) {
					bestBaseParams[i] = baseParams[i];
				}
			
				for (int i = 0; i < pownum; i++)
				for (int j = 0; j<AnnotAlphabetSize; j++)
				{ 
					bestAnnotParams[i][j] = annotParams[i][j];  
				}	  
			}
			
		}   // closes the for loop that tries different initializations
		/************SORT THE PARAMS*****************************************************************/

		for (int i=0; i<pownum; i++) {
			baseParamsPrev[i] = bestBaseParams[i];
		}
		for (int i = 0; i < pownum; i++)
		for (int j=0; j<AnnotAlphabetSize; j++) {
		  annotParamsPrev[i][j] = bestAnnotParams[i][j];
		}
		parAPrev = bestParA;
		parBPrev = bestParB;
		biasSPrev = bestBiasS;
		biasAPrev = bestBiasA;

		/****************  FPOUT *****************************************************************/  
		outParam<<"Motif width "<<mw<<endl;
		outParam<<"Base parameters "<<endl;
		//scale base parameters so that a logo can be generated
		double sum;
		double maxValue = -1000;	
		for (int i=0; i<pownum; i++) {
			if(bestBaseParams[i] > maxValue)
				maxValue = bestBaseParams[i];
		}
		double cons = (-1 * biasS) / maxValue;  // so that the best motif will have affinity 0.5
						
		for (int i=0; i<pownum; i++) {
			bestBaseParamsScaled[i] = bestBaseParams[i] * cons;
		}

		// print the scaled parameters
		for (int i=0; i<pownum; i++) {
			outParam<<getString(i, motifWidth)<< "\t" <<-1 * bestBaseParamsScaled[i];
			for (int j = 0; j < AnnotAlphabetSize; j++)
				outParam << "\t" << bestAnnotParams[i][j];
			outParam << "\n";
		}

		double* relativeAffinities =  (double *)malloc(sizeof(double) * AnnotAlphabetSize);
		double maxAffinity = 0;
		int index = 0;
		for (int i=0; i<AnnotAlphabetSize; i++) {
			relativeAffinities[i] = 1 / (1 + exp(-1 * (mw * bestAnnotParams[index][i] + biasA)));
			if( relativeAffinities[i] > maxAffinity )
				maxAffinity = relativeAffinities[i];
		}
		outParam<<endl<<"Relative affinities to the structural contexts "<<endl;
		for (int i=0; i<AnnotAlphabetSize; i++) {
			relativeAffinities[i] = relativeAffinities[i] / maxAffinity;
			outParam<<annotAlphabet[i]<<"\t"<<relativeAffinities[i]<<endl;
		}
		free(relativeAffinities);
		outParam<<endl;

				
		/****************  PRINTING PARAMETERS *****************************************************************/  
		//assign the best parameters to the current parameters so that they are used in the best kmers calculation
		
		for (int i=0; i<pownum; i++) {
			baseParams[i] = bestBaseParams[i];
		}
		
		for (int i = 0; i < pownum; i++)
			for (int j = 0; j < AnnotAlphabetSize; j++) {
				annotParams[i][j] = bestAnnotParams[i][j];
			}
											
		parA = bestParA;
		parB = bestParB;
		biasS = bestBiasS;
		biasA = bestBiasA;
		
		//calculate the scores with the best parameters
		for(int i=0;i<trainSeqNum;i++)
		{
			scoresTrain[i] = calcScore(i);			
		}
		
		
		
		if(AnnotAlphabetSize > 1) //alph 1 2 or 3 
			numParams = pownum - 1 + AnnotAlphabetSize - 1 + 2 + 2;
		else 
			numParams = pownum - 1 + AnnotAlphabetSize + 2 + 2;
		
		cout << "FILES" << endl;
		
		modelFileName = dir + "/model_" + outFile + "_" + mwStr + ".txt" ;
		testFileName = dir +"/test_" + outFile + "_" +  mwStr + ".txt" ;
		testFileNameSeq = dir + "/test_" + outFile + "_" +  mwStr + "_seq.txt" ;
		rankFileName = dir + "/train_" + outFile + "_" +  mwStr + ".txt" ;

		/* print the training sequences with their ratios and scores to a file */
		outTrain.open(rankFileName.c_str());
		for (int i=0; i < trainSeqNum; i++)
		{
   		  	outTrain<<ratiosTrain[i]<<"\t"<<scoresTrain[i]<<endl;  	
		}	
		outTrain.close();
		
		/* use the best obtained parameters to calculate the error on the tests set */
		testBestParams(testFileName, modelFileName, numParams); 
		testBestParamsSeq(testFileNameSeq, modelFileName, numParams); 

		modelOut.open(modelFileName.c_str(), ios::app);
		modelOut<<"Error on the training set: "<<bestError<<endl;
		modelOut<<"Number of parameters: "<<numParams<<endl<<endl;
		modelOut<<"Base Parameters"<<endl;
		
		for (int i=0; i<pownum; i++) {
			modelOut<<bestBaseParams[i];
			for (int j = 0; j < AnnotAlphabetSize; j++)
				modelOut << "\t" << bestAnnotParams[i][j];
			modelOut<<"\n";
		}
		
		modelOut<<endl;
		modelOut<<"Annot Parameters for each context in the alphabet: "<<annotAlphabet<<endl;
		for (int i=0; i<AnnotAlphabetSize; i++) {
			modelOut<<annotAlphabet[i]<<"\t"<<bestAnnotParams[index][i] <<endl;
		}
		modelOut<<endl<<endl;
		modelOut<<"Beta_s (bias in sequence model) : "<<bestBiasS<<endl;
		modelOut<<"Beta_p (bias in structure context model) : " <<bestBiasA<<endl;
		modelOut<<"alpha (scaling factor) : "<<bestParA<<endl;
		modelOut<<"b (bias in least squares error model) : "<<bestParB<<endl;
		modelOut.close();
		// TO DO				
		findBestScores(modelFileName, motifWidth, 20); //finds the top 20 kmers and align their profile
		// to get position specific structure preferences
		// free everything
		free(baseParams);
		free(bestBaseParams);
		free(bestBaseParamsScaled);
		// delta values
		for (int i = 0; i < trainSeqNum; i++)
		{
			for (int k = 0; k < AnnotAlphabetSize; k++)
				free(sumAnnotParamsTrain[i][k]);
			free(sumAnnotParamsTrain[i]);

		}
		
	}   // CLOSES THE BIG LOOP --MOTIFWIDTH CHANGES
	//free  everything that is allocated in main
	for(int i = 0 ; i < trainSeqNum ; i++)
	{
		free(sumBaseParamsTrain[i]);
		free(AffinitySeqVals[i]);	
		free(AffinityAnnotVals[i]);
		free(affinityValsTrain[i]);	
	}
	free(baseParamsPrev);

	free(means);
	free(variances);
	free(sumBaseParamsTrain);
	free(sumAnnotParamsTrain);
	for (int i = 0; i < pow(AlphabetSize, maxWidth); i++) free(annotParams[i]);
	free(annotParams);
	free(AffinitySeqVals);
	free(AffinityAnnotVals);
	free(affinityValsTrain);
	free(kmerNumsTrain);
	free(kmerNumsTest);
}
// Once the test sequences are scored, find the best 20 kmer scores and output the average annotation profile of these 20 kmers
void findBestScores(string filename, int motifWidth, int numBest)
{
	ofstream out;
  	out.open(filename.c_str(), ios::app);
	int count = 0;
	int kmerCount = seqsTrain[0].length()- motifWidth +1 ;
	
	vector< vector <double> > profiles(AnnotAlphabetSize); 
	for(int l = 0 ; l< AnnotAlphabetSize ; l++)
    {
	 	profiles[l].resize(motifWidth, 0);
    }
	
	vector< kmerScores> allKmerScores; //(trainSeqNum * kmerCount);

	for(int seqID = 0; seqID < trainSeqNum ; seqID++)
	{
		kmerCount = seqsTrain[seqID].length() - motifWidth + 1;
		for(int kmer= 0 ; kmer < kmerCount ; kmer++)	
		{
			kmerScores kscore;
			kscore.seqID = seqID;
			kscore.kmerNum = kmer;
			kscore.affinity = affinityValsTrain[seqID][kmer];
			kscore.seqAffinity = AffinitySeqVals[seqID][kmer];
			kscore.annotAffinity = AffinityAnnotVals[seqID][kmer];
			allKmerScores.push_back (kscore);
			count++;
		}
    }
	
    sort(allKmerScores.begin(), allKmerScores.end(), sort_by_one());
	out<<endl;
	out<<"Top 20 kmers: \t sequence number \t kmer number"<<endl;
    for(int best = 0 ; best < numBest; best++)
    {
	    out<<setw(12);
    	for(int m = 0 ; m < motifWidth ; m++)
    	{	
			out<<seqsTrain[allKmerScores[best].seqID][allKmerScores[best].kmerNum+m];
     	}
	    out<<setw(15)<<allKmerScores[best].seqID<<setw(15)<<allKmerScores[best].kmerNum<<endl;   	
    }
	out<<endl;
	out<<"Averaged annotation profiles of the top 20 kmers"<<endl;
    for(int l = 0 ; l< AnnotAlphabetSize ; l++)
    {  	 
    	out<<annotAlphabet[l]<<"\t"; 	    	 	 
     	for(int m = 0; m < motifWidth; m++)
     	{
			for(int best = 0 ; best < numBest; best++)
			{
				profiles[l][m] += inputAnnotsTrain[allKmerScores[best].seqID][l][allKmerScores[best].kmerNum+m];
			}
			profiles[l][m] /= numBest;
			out<<profiles[l][m]<<"\t";			
	    }		
        out<<endl;
    } 
    out.close();
    allKmerScores.clear();
    profiles.clear();    
}



// Calculates the score of the sequence with sequence id seqID
double calcScore(int seqID)
{
  double score = 0;
  double storeSumBase = 0;
  double sumAnnotWeights = 0 ;
  
  for(int kmer = 0; kmer< kmerNumsTrain[seqID]; kmer++) // all kmers
  {
  		sumBaseParamsTrain[seqID][kmer] = sumBaseParamsFunc(seqID, kmer); 		
  }
  for(int kmer = 0; kmer< kmerNumsTrain[seqID] ; kmer++) // all kmers
  {
  	 	storeSumBase = sumBaseParamsTrain[seqID][kmer];
  	 	AffinitySeqVals[seqID][kmer] = 1.0 / (1.0 + exp( -1 * (biasS + storeSumBase))); 
  	 	sumAnnotWeights = 0;
  	 	for(int letter = 0; letter < AnnotAlphabetSize ; letter++) // for L,R,P,U,H etc.
  	 	{
  	 	    sumAnnotParamsTrain[seqID][letter][kmer] = sumProfilesFunc(seqID, letter, kmer);
			int index = getInt(seqsTrain[seqID].c_str(), kmer, motifWidth) ;
  	 	    sumAnnotWeights += annotParams[index][letter] * sumAnnotParamsTrain[seqID][letter][kmer];   
	        }
		AffinityAnnotVals[seqID][kmer] = 1.0 / (1.0 + exp(-1 * (biasA + sumAnnotWeights)));
	    affinityValsTrain[seqID][kmer] = (AffinitySeqVals[seqID][kmer] * AffinityAnnotVals[seqID][kmer]) - SMALL1;
   	 	
		// note the score here is simply a sum of all binding probabilities of each k-mer, which is 
		//  different from the RCK paper where score=1-prod_{all k-mers of the seq}(1-affinityValsTrain[seqID][kmer]). 
		//  Also, note RNAcontext's implementation is also like this which is different from RNAcontext paper
   	 	score +=  affinityValsTrain[seqID][kmer];		
}
  
  return score; // / kmerNumsTrain[seqID];	// change to be same as fix2
}

double calcScoreTestWindow(int seqID, int window)
{
  double score = 0;
  double storeSumBase = 0;
  double sumAnnotWeights = 0 ;
  double sumAnnot;
  double affSeq, affAnnot, occ;
  double max = 0;

  for(int kmer = 0; kmer< kmerNumsTest[seqID] ; kmer++) // all kmers //
  {
                storeSumBase = sumBaseParamsTestFunc(seqID, kmer);
                affSeq = 1.0 / (1.0 + exp( -1 * (biasS + storeSumBase))); 
                sumAnnotWeights = 0; //reset sumAnnotWeights 
                for(int letter = 0; letter < AnnotAlphabetSize ; letter++) // for L,R,P,U,H etc.
                {
                         int index = getInt(seqsTest[seqID].c_str(), kmer, motifWidth) ;
                        sumAnnot = sumProfilesFuncTest(seqID, letter, kmer);
                    sumAnnotWeights += annotParams[index][letter] * sumAnnot;
            }
            affAnnot = 1.0 / (1.0 + exp(-1 * (biasA + sumAnnotWeights)));
            occ = (affSeq * affAnnot);// - small1;
            score += occ;
	if (occ > max) max = occ;
  }

  return score; // / (kmerNumsTest[seqID]);	// change to be same as fix2
}

// Calculates the score of a test sequence
double calcScoreTest(int seqID)
{
  double score = 0;
  double storeSumBase = 0;
  double sumAnnotWeights = 0 ;
  double sumAnnot;
  double affSeq, affAnnot, occ;

  for(int kmer = 0; kmer< kmerNumsTest[seqID] ; kmer++) // all kmers //
  {
  	 	storeSumBase = sumBaseParamsTestFunc(seqID, kmer);
  	 	affSeq = 1.0 / (1.0 + exp( -1 * (biasS + storeSumBase))); 
  	 	sumAnnotWeights = 0; //reset sumAnnotWeights 
  	 	for(int letter = 0; letter < AnnotAlphabetSize ; letter++) // for L,R,P,U,H etc.
  	 	{	
			 int index = getInt(seqsTest[seqID].c_str(), kmer, motifWidth) ;
  	 		sumAnnot = sumProfilesFuncTest(seqID, letter, kmer);
  	 	    sumAnnotWeights += annotParams[index][letter] * sumAnnot;   
	    }
	    affAnnot = 1.0 / (1.0 + exp(-1 * (biasA + sumAnnotWeights)));
	    occ = (affSeq * affAnnot);// - small1;

		// note the score here is simply a sum of all binding probabilities of each k-mer, which is 
		//  different from the RCK paper where score=1-prod_{all k-mers of the seq}(1-affinityValsTrain[seqID][kmer]). 
		//  Also, note RNAcontext's implementation is also like this which is different from RNAcontext paper
		score += occ;
  }

  return score; // / (kmerNumsTest[seqID]);	// change to be same as fix2
}
// Calculates the score of a test sequence
double calcScoreTestSeq(int seqID)
{
  double score = 0;
  double storeSumBase = 0;
  double affSeq, occ;

  for(int kmer = 0; kmer< kmerNumsTest[seqID] ; kmer++) // all kmers //
  {
  	 	storeSumBase = sumBaseParamsTestFunc(seqID, kmer);
  	 	affSeq = 1.0 / (1.0 + exp( -1 * (biasS + storeSumBase))); 
		occ = affSeq;// - small1;
   	 	score += occ;
  } 

  return score;
}

// Calculate the gradients
double calcGradParams(int paramType, int alphIndex = 0, ap::real_1d_array& gradVec = AVAL)
{
  double kmerVal = 0;
  double storeResult = 0;
  double storeMult = 0;
  double result = 0; 
  double initialTerm = 0;
  int deltaVal;

  if (paramType == BASE) {
    for (int i = 0; i < alphIndex; i++) {
      gradVec(i+1) = 0;
    }
  }

  if (paramType == ANNOT) {
    for (int i = 0; i < alphIndex; i++) 
	for (int j = 0; j < AnnotAlphabetSize; j++)
		gradVec(i*AnnotAlphabetSize + j + 1) = 0;// comment out: //gradVec(i*AnnotAlphabetSize+j+alphIndex+1) = 0;
  }

  for(int seqID	= 0; seqID < trainSeqNum; seqID++)
  { 		
		initialTerm = 2 * (ratiosTrain[seqID] - ( parA * scoresTrain[seqID]) - parB); // 2 * sum_i (x_i - a*. f_i - b* )
		
		if(paramType == A)
		{
			storeResult += initialTerm * -1 * scoresTrain[seqID];
		}
		else if(paramType == B)
		{
			storeResult += initialTerm * -1;
		}
		else
		{
			for(int kmer = 0; kmer< kmerNumsTrain[seqID] ; kmer++) // for each kmer
			{
  	 			if (paramType == BASE) {
					kmerVal = (AffinityAnnotVals[seqID][kmer] * AffinitySeqVals[seqID][kmer] * (1-AffinitySeqVals[seqID][kmer]));
					int index = getInt(seqsTrain[seqID].c_str(), kmer, motifWidth);
					gradVec(index+1) = gradVec(index+1) + kmerVal * initialTerm * -1 * parA;
				}
				else if (paramType == ANNOT)
				{
					int index = getInt(seqsTrain[seqID].c_str(), kmer, motifWidth);
					for (int j = 0; j < AnnotAlphabetSize; j++) {
						kmerVal = (AffinitySeqVals[seqID][kmer] * AffinityAnnotVals[seqID][kmer] * (1 - AffinityAnnotVals[seqID][kmer]) * sumAnnotParamsTrain[seqID][j][kmer]);
						gradVec(index*AnnotAlphabetSize + j + 1) = gradVec(index*AnnotAlphabetSize + j + 1) + kmerVal * initialTerm * -1 * parA;// comment out: //gradVec(alphIndex + index*AnnotAlphabetSize + j + 1) = gradVec(alphIndex + index*AnnotAlphabetSize + j + 1) + kmerVal * initialTerm * -1 * parA;
					}
				}
        	    else if(paramType == BIASS) 
        	    	 kmerVal +=  AffinityAnnotVals[seqID][kmer] * AffinitySeqVals[seqID][kmer] * (1-AffinitySeqVals[seqID][kmer]);    
        	    else if(paramType == BIASA) 
        	     	 kmerVal += (AffinitySeqVals[seqID][kmer] * AffinityAnnotVals[seqID][kmer] * (1 - AffinityAnnotVals[seqID][kmer]));           	
            }
            storeResult += initialTerm * -1 * parA * kmerVal ;  // negative log likelihood
            kmerVal = 0;
       }
  }

   result = storeResult ;  
   if(paramType == A)
		result = parA * result;
   if (paramType == BASE) {
	   for (int i = 0; i < alphIndex; i++) {
		   double tmp = exp(-biasS - baseParams[i]);
	   }
   }
   if (paramType == BIASS) {		// this is the gradient for BiasS in the regularization term of RCK
	   for (int i = 0; i < alphIndex; i++) {
		   double tmp = exp(-biasS - baseParams[i]);
		   result += SMALL2 / pow(1 + tmp, 2) * tmp;
	   }
   }
   else if (paramType == ANNOT)
	   for (int i = 0; i < alphIndex; i++)
		   for (int j = 0; j < AnnotAlphabetSize; j++) {
		   }

   return result;
}
//calculate the gradients and store them in gradVec
void gradientFunc(ap::real_1d_array& gradVec, ap::real_1d_array& paraVec) 
{
    int count = 1;
	// map back the base parameters
	int pownum = pow(AlphabetSize, motifWidth);
	// comment out: //calcGradParams(BASE, pownum, gradVec);   count += pownum;

   // map back the annot parameters
	calcGradParams (ANNOT,pownum, gradVec);	  
	count+= pownum*AnnotAlphabetSize;
	//map back bias
	gradVec(count) = calcGradParams(BIASS, pownum); 	
	count++;
	
	gradVec(count) = calcGradParams(BIASA, pownum);  
	count++;
	//map back A
	gradVec(count) = calcGradParams(A);	
	count++;
	//map back B
	gradVec(count) = calcGradParams(B);		
}
void funcgrad( ap::real_1d_array& x, double& f, ap::real_1d_array& g)
{
	f=0;
	f = originalFunc(x);
	gradientFunc(g, x);		// note gradient vector g has the same dimension as parameter vector x
	iter++;
}


// No motif model is learned in this case, given a motif, test sequences are scored
void scoreSeqs(char *dataFileName, char *annotFileName, char *motifFileName, char *outDir, int width, int ignoreAnnot)
{
	string dir = outDir;
	string motifFileNameStr = motifFileName;	
                string mwStr;
                stringstream outs;
                outs << width;
                mwStr = outs.str();
	int txtIndex = motifFileNameStr.find(".txt");
	int lastSlash = motifFileNameStr.rfind("/") + 1;
	//CHANGE THIS
	testSeqNum = readData(dataFileName,seqsTest, ratiosTest);
	kmerNumsTest = (int *) malloc(sizeof(int) * testSeqNum);
	
	vector<double> stds(AnnotAlphabetSize);
	// Read the motif model from file
	ifstream inputMotif;
	string modelFileName = dir + "/model_" + motifFileNameStr + "_" + mwStr + ".txt"; 
	inputMotif.open(modelFileName.c_str());
	cout << "motifFileName = " << modelFileName << endl;
	if(inputMotif.fail())
		cout<<" cannot open motif filename"<<endl;

	ofstream out;
	inputAnnotsTest.resize(testSeqNum);
	
	motifWidth = width;
	for(int i = 0 ; i < testSeqNum ; i++)
	{	
		kmerNumsTest[i] = seqsTest[i].length() - motifWidth +1;
	}
	
    for(int i = 0 ; i< testSeqNum ; i++)
    {
		inputAnnotsTest[i].resize(AnnotAlphabetSize);	 

		for(int k = 0; k <AnnotAlphabetSize ; k++)
			inputAnnotsTest[i][k].resize(seqsTest[i].length()); 
    }
	readAnnots(annotFileName, inputAnnotsTest, testSeqNum, ignoreAnnot, false); //testSeqNum is assigned here
	//score the sequences
    baseParams = (double *)malloc(pow(AlphabetSize, motifWidth)*sizeof(double));
    annotParams = (double **)malloc(pow(AlphabetSize, motifWidth)*sizeof(double*));
	for (int i = 0; i< pow(AlphabetSize, motifWidth); i++)
		annotParams[i] = (double *)malloc(AnnotAlphabetSize*sizeof(double));
 
	string line, foo;	
	double dataTemp;

	// read the base parameters --Theta
	for(int i = 0; i < 5; i++) getline(inputMotif, line);
	int pownum=pow(AlphabetSize, motifWidth);
	for(int i = 0; i < pow(AlphabetSize, motifWidth); i++)
	{
		getline(inputMotif, line);
		istringstream ins;
		ins.str(line);
			ins>>baseParams[i];
		for (int j=0; j<AnnotAlphabetSize; j++) 
		{
			ins>>dataTemp;
			annotParams[i][j] = dataTemp;
		}
	}
	for (int i = 0; i < 4 + AnnotAlphabetSize; i++) getline(inputMotif,line);

	istringstream ins; string delimiter = ":";
        getline(inputMotif, line); ins.str(line.substr(line.find(delimiter)+2,10));
	ins >> biasS; ins.clear();
        getline(inputMotif, line); ins.str(line.substr(line.find(delimiter)+2,10));
        ins >> biasA; ins.clear();
        getline(inputMotif, line); ins.str(line.substr(line.find(delimiter)+2,10));
        ins >> parA; ins.clear();
        getline(inputMotif, line); ins.str(line.substr(line.find(delimiter)+2,10));
        ins >> parB; ins.clear();
	cout << biasS << " " << biasA << " " << parA << " " << parB << endl;

        vector< kmerScores> allKmerScores; //(trainSeqNum * kmerCount);
        cout << "Scoring k-mers" << endl;       
        for(int i = 0; i < pownum; i++)
        {
                for(int a = 0 ; a < AnnotAlphabetSize ; a++)    
                {
                        kmerScores kscore;
                        kscore.seqID = i;
                        kscore.kmerNum = a;
                        kscore.affinity = 1.0 / ( (1.0 + exp( -1 * (biasS + baseParams[i]))) * (1.0 + exp( -1 * (biasA + annotParams[i][a]))) );
	                kscore.seqAffinity = annotParams[i][a];
                        kscore.annotAffinity = baseParams[i];
                        allKmerScores.push_back (kscore);
                }
	    }
	    sort(allKmerScores.begin(), allKmerScores.end(), sort_by_one());
	double* params = (double*)malloc(sizeof(double) * AnnotAlphabetSize);
	double** pwm = (double**)malloc(sizeof(double*) * width);
	for (int i =0; i < width; i++)
		pwm[i] = (double*)malloc(sizeof(double) * 4);

	int max = allKmerScores[0].kmerNum;
	double sum=0;
	for (int i = 0; i < pownum * AnnotAlphabetSize; i++)
		if (allKmerScores[i].seqID == allKmerScores[0].seqID) {
			params[allKmerScores[i].kmerNum] = allKmerScores[i].affinity;
			sum += params[allKmerScores[i].kmerNum];
		}

	for (int i = 0; i < AnnotAlphabetSize; i++)
		cout << annotAlphabet[i] << ":\t" << params[i] / sum << endl;

	string con = getString(allKmerScores[0].seqID, width);
	for (int i = 0; i < pownum * AnnotAlphabetSize; i++) {
		int match = 0;
		string tmp = getString(allKmerScores[i].seqID, width);
		int dif = -1;
		for (int t = 0; t < width; t++)
			if (tmp[t] == con[t]) match++;
			else dif = t;
		if (match == width - 1 && allKmerScores[i].kmerNum == max) {
			pwm[dif][basetoInt(tmp[dif])] = allKmerScores[i].affinity;
		}
	}
	for (int i = 0; i < width; i++) {
		pwm[i][basetoInt(con[i])] = allKmerScores[0].affinity;
		double sum = 0;
		for (int a = 0; a < 4; a++) sum += pwm[i][a];
		for (int a = 0; a < 4; a++) pwm[i][a] /= sum;
	}
	for (int a = 0; a < 4; a++) {
		cout << getString(a,1) + ":";
		for (int i = 0; i < width; i++)
			cout << "\t" << pwm[i][a];
		cout << endl;
	}

	
               for(int i = 0 ; i < testSeqNum ; i++)
                {       
                        kmerNumsTest[i] = seqsTest[i].length()- motifWidth +1; 
	}
	scoresTest.resize(testSeqNum);
        for(int i=0;i<testSeqNum;i++) {
            scoresTest[i] = calcScoreTestWindow(i, 1);
	}
		string rankFileName = dir + "/pred_" + motifFileNameStr + "_" + mwStr + ".txt"; 
		ofstream outTest;
                outTest.open(rankFileName.c_str());
                for (int i=0; i < testSeqNum; i++)
                {
                        outTest<<ratiosTest[i]<<"\t"<<scoresTest[i]<<endl;   
                }       
                outTest.close();
}

int readModel(string inputFileName, double *baseParams, double **annotParams, double &biasS, double &biasA, double &alpha, double &beta, int pownum, int size)
{
	int counter = 0;
        istringstream ins;
        double ratio;
        string seq, line, foo;
        ifstream input;
        input.open(inputFileName.c_str());
        if(input.fail()) {
               cout<<" cannot open the input file "<<endl;
        }
        else
        {
		// read header
		for (int i = 0; i < 5; i++) getline(input, line);

		for (int i = 0; i < pownum; i++) {
			getline(input, line);
			istringstream ins;
			ins.str(line);
			ins >> baseParams[i];
			for (int j = 0; j < size; j++)
				ins >> annotParams[i][j];
		}

		for (int i = 0; i < 8; i++) getline(input, line);
		getline(input, line, ':'); getline(input, line, ':'); ins.str(line);
		ins >> biasS;
                getline(input, line, ':'); ins.str(line);
                ins >> biasA;
                getline(input, line, ':'); ins.str(line);
                ins >> alpha;
                getline(input, line, ':'); ins.str(line);
                ins >> beta;
		counter++;		
        }
        return counter;

}

// Reads the input file, fills the global variables seqs and ratios, returns the numeber of probes in the file.
int readData(string inputFileName, vector<string> & seqs, vector<double> & ratios)
{
	int counter = 0;
	istringstream ins;
	double ratio;
	string seq, line, foo;
	ifstream input;
	input.open(inputFileName.c_str());
	if(input.fail())
	{
		cout<<" cannot open the input file "<<endl;
	}
	else
	{
		while(getline(input,line))
		{
			istringstream ins;
			ins.str(line);
			ins>>ratio>>seq;
			replace(seq.begin(), seq.end(), 'T','U');
			seqs.push_back(seq);
			ratios.push_back(ratio);
			counter++;
		}
	}
	return counter;
}	

/*       MAIN        */
int main(int argc, char **argv)
{
	
	InputParams pm;
	// set the parameters
	setParameters(argc, argv, &pm);   // annotalphabetsize is determined here
	// open the log file
	openlogfile(pm.logFileName, MAX_LOGENTRIES);
	
	cout << "open log file" << endl;
	// convert the alphabet to uppercase 
	Upper(pm.alphabet);  
	if (pm.motifFileNameDefinedFlag == 0)	// i.e. training mode
	{

		trainSeqNum = readData(pm.trainDataFileName, seqsTrain, ratiosTrain); //read the file
		testSeqNum = readData(pm.testDataFileName, seqsTest, ratiosTest); //read the file		

		// Memory allocation for global variables.
		scoresTrain.resize(trainSeqNum);
		scoresTest.resize(testSeqNum);
		kmerNumsTrain = (int *)malloc(sizeof(int) * trainSeqNum);
		kmerNumsTest = (int *)malloc(sizeof(int) * testSeqNum);
		affinityValsTrain = (double **)malloc(sizeof(double*) * trainSeqNum);
		AffinitySeqVals = (double **)malloc(sizeof(double*) * trainSeqNum);
		AffinityAnnotVals = (double **)malloc(sizeof(double*) * trainSeqNum);
		sumBaseParamsTrain = (double **)malloc(sizeof(double*) * trainSeqNum);
		sumAnnotParamsTrain = (double ***)malloc(sizeof(double**) * trainSeqNum);
		means = (double *)malloc(sizeof(double) * AnnotAlphabetSize);
		variances = (double *)malloc(sizeof(double) * AnnotAlphabetSize);

		inputAnnotsTrain.resize(trainSeqNum);
		inputAnnotsTest.resize(testSeqNum);
		for (int i = 0; i < trainSeqNum; i++)
		{
			scoresTrain[i] = 0;
			kmerNumsTrain[i] = 0;
			inputAnnotsTrain[i].resize(AnnotAlphabetSize);
			for (int k = 0; k < AnnotAlphabetSize; k++)
				inputAnnotsTrain[i][k].resize(seqsTrain[i].length());

		}
		for (int i = 0; i < testSeqNum; i++)
		{
			inputAnnotsTest[i].resize(AnnotAlphabetSize);
			for (int k = 0; k < AnnotAlphabetSize; k++)
				inputAnnotsTest[i][k].resize(seqsTest[i].length());

		}
		// Read the profiles for training and test sequences
		readAnnots(pm.annotFileNameTrain, inputAnnotsTrain, trainSeqNum, pm.ignoreAnnotFlag, true);
		readAnnots(pm.annotFileNameTest, inputAnnotsTest, testSeqNum, pm.ignoreAnnotFlag);
		// Start searching motifs	// note searchMotif() has one more parameter now for external PWM file name
		searchMotif(pm.seedNum, pm.seedWidth, pm.minWidth, pm.maxWidth, pm.outFileName, pm.annotFileNameTrain, pm.outDir, pm.maxIter, pm.inputModelFlag, pm.pwmFileName);
	}
	else // score the test sequences given a motif model 	// i.e. predicting mode
	{
		scoreSeqs(pm.testDataFileName, pm.annotFileNameTest, pm.motifFileName, pm.outDir, pm.minWidth, pm.ignoreAnnotFlag);
	}
	
	return 0;
}

