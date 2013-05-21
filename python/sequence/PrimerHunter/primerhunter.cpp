using namespace std;

#include "primerhunter.h"
#include "stdafx.h"
#include <string>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "nnparams.h"
#include "thermalign.h"
#include "dinkelbach.h"
#include "dpal.h"

#include <stdio.h>

char * targets [TARGETS_MAX_SIZE];
int nTargets=0;
int seqPositions=0;
int coveredTargets[TARGETS_MAX_SIZE]={0};
int nCoveredTargets = 0;
char * nontargets [NON_TARGETS_MAX_SIZE];
int nNonTargets=0;
TSeed seedsTableTargets [SEEDSTABLE_TAM];
int nSeedTargets;
TSeed seedsTableNonTargets [SEEDSTABLE_TAM];
int nSeedNonTargets;
Primer primers [CANDIDATES_MAX_SIZE];
int nPrimers=0;
Pair pairs [CANDIDATES_MAX_SIZE];
int nPairs = 0;
//Place where primers loaded from primers file are stored
Primer loadedPrimers[CANDIDATES_MAX_SIZE];
int nLoadedPrimers=0;


//User parameters
char targetsFile [1000]="";
char nonTargetsFile [1000]="";
char primersFile [1000]="";
char primersOutputFile [1000]="";
int minLength = DEF_MIN_PRIMER_LENGTH;
int maxLength = DEF_MAX_PRIMER_LENGTH;
int productMinLength = DEF_MIN_PROD_LENGTH;
int productMaxLength = DEF_MAX_PROD_LENGTH;
char forwardPrimer [PRIMERS_MAX_SIZE] = DEF_FORWARD_PRIMER;
char reversePrimer [PRIMERS_MAX_SIZE] = DEF_REVERSE_PRIMER;
int beginPosForward = DEF_BEGINPOS_FORWARD;
int endPosForward = DEF_ENDPOS_FORWARD;
int beginPosReverse = DEF_BEGINPOS_REVERSE;
int endPosReverse = DEF_ENDPOS_REVERSE;
char maskTargets [PRIMERS_MAX_SIZE] = DEF_MASK_TARGETS;
char maskNonTargets [PRIMERS_MAX_SIZE] = DEF_MASK_NONTARGETS;
char degeneracy[PRIMERS_MAX_SIZE]=DEF_DEGENERACY;
int nseqForCandidates = DEF_SEQ_FOR_CANDS;
float minCoverageTargets = DEF_MIN_COVERAGE_TARGETS;
float maxCoverageNonTargets = DEF_MAX_COVERAGE_NONTARGETS;
float maxSelfScore = DEF_MAX_SELF_SCORE;
float maxEndScore = DEF_MAX_END_SCORE;
float minGCContent = DEF_MIN_GCCONTENT; 
float maxGCContent = DEF_MAX_GCCONTENT;
int gcClamp = DEF_GCCLAMP;
int maxPolyX = DEF_MAX_POLY_X;
float concPrimers = DEF_CONC_PRIMERS; 
float concSequences = DEF_CONC_SEQUENCES; 
float salt = DEF_SALT;
int salMethod = SALT_METHOD_SANTALUCIA; 
float minTempTargets=DEF_MIN_TEMP_TARGETS;
float maxTempTargets=DEF_MAX_TEMP_TARGETS;
float maxTempNonTargets=DEF_MAX_TEMP_NONTARGETS;
float deltaTempNonTargets=DEF_DELTA_TEMP_NONTARGETS;
float maxPairTempDiff=DEF_MAX_PAIR_TEMP_DIFF;
char primersLabel [1000] = DEF_PRIMERS_LABEL;
int fullStats = DEF_FULL_STATS;
char seqRank [10000] = DEF_SEQ_RANK;

//Melting temperature variables
PNNParams paramsPrimerSeq;

// dpal variables

dpal_args local_args;
dpal_args end_args;

//General methods over sequences

int getSequenceMaxLength (char * sequences [], int nseq) {
	int answer = 0;
	int max = 0;
	for (int i=0;i<nseq;i++) {
		int l = strlen(sequences[i]);
		if(l>max) {
			max = l;
			answer = i;
		}
	}
	return answer;
}

int getHash (char * seed) {
	int number =0;
	for(int i=0;seed[i]!=0;i++) {
		number*=4;
		if(seed[i]=='C') {
			number+=1;
		}
		if(seed[i]=='G') {
			number+=2;
		}
		if(seed[i]=='T') {
			number+=3;
		}
	}
	return number;
}

void buildSeed (char * answer, int number,int size) {
	answer[size] = 0;
	for(int i=0;i<size;i++) {
		int nextDigit = number%4;
		int index = size-i-1;
		if(nextDigit == 3){
			answer[index]='T';
		} else if (nextDigit == 2) {
			answer[index]='G';
		} else if (nextDigit == 1) {
			answer[index]='C';
		} else {
			answer[index]='A';
		}
		number = number/4;
	}
}

void buildSeedFromSequence (char * seq, char * mask, char * seed) {
	int lmask = strlen(mask);
	int lseq = strlen (seq);
	int seedIndex = 0;
	for(int i=0; i < lmask && i < lseq; i++) {
		if(mask[i] == '1') {
			seed[seedIndex] = seq[lseq-i-1];
			seedIndex++;
		}
	}
	seed[seedIndex] = 0;
}

char getComplement(char base) {
	if(base == 'A') return 'T';
	else if(base == 'G') return 'C';
	else if(base == 'C') return 'G';
	else if(base == 'T') return 'A';
	else if(base == 'N') return 'N';
	return 'A';
}

float countGCContent(char * seq ) {
	int lseq = strlen(seq);
	float count = 0;
	//Check degeneracy
	for(int k=0;k<lseq;k++) {
		if (seq[k] == 'N') {
			count+=0.5;
		}
		if (seq[k] == 'G' || seq[k] == 'C' ) {
			count+=1;
		}
	}
	return count;
	
}

void calculateComplement(char * answer, char * sequence) {
	int nseq = strlen(sequence);
	for (int i= 0; i<nseq; i++) {
		answer[i] = getComplement(sequence[i]);
	}
	answer [nseq] = 0;
}

void calculateReverse(char * answer, char * sequence) {
	int nseq = strlen(sequence);
	for (int i= nseq-1; i>=0; i--) {
		answer[nseq-i-1] = sequence[i];
	}
	answer [nseq] = 0;
}

void calculateReverseComplement(char * answer, char * sequence) {
	int nseq = strlen(sequence);
	for (int i= nseq-1; i>=0; i--) {
		answer[nseq-i-1] = getComplement(sequence[i]);
	}
	answer [nseq] = 0;
}

void reverseComplementSequences (char * sequences [], int n) {
	char next [SEQUENCES_MAX_SIZE];
	for(int i=0; i< n; i++ ) {
		calculateReverseComplement(next, sequences[i]);
		strcpy(sequences[i],next);
	}
}
void addDollarSigns(char * seq) {
	int n = strlen(seq);
	seq[n+2]=0;
	seq[n+1]='$';
	for(int i=n;i>0;i--) {
		seq[i]=seq[i-1];
	}
	seq[0] = '$';
}
void removeDollarSigns(char * seq) {
	int n = strlen(seq);
	for(int i=0;i<n-2;i++) {
		seq[i]=seq[i+1];
	}
	seq[n-2]=0;
}
/**
 * seq1 5' - 3'
 * seq2 3' - 5'
 */
float calculateMeltingTemperature (char * seq1, char * seq2, PNNParams params) {
	int length = strlen(seq2);
	int perfectComplements = 1;
	if(strlen(seq1)!= length) {
		perfectComplements = 0;
	} 
	for(int k=0;k<length && perfectComplements;k++) {
		if(seq2[k] != getComplement(seq1[k])) {
			perfectComplements = 0;
		}
	}
	addDollarSigns(seq1);
	addDollarSigns(seq2);
	CThermAlign myAlignment(strlen(seq1),strlen(seq2),params);
	myAlignment.InitStrings(seq1,seq2,strlen(seq1),strlen(seq2));
	myAlignment.InitBorder();
	myAlignment.CalculateTable();
	float answer;
	if(perfectComplements) {
		answer = myAlignment.GetMeltingTempC(myAlignment.maxloci,myAlignment.maxlocj);
	} else {
		float tempK = myAlignment.GetMeltingTempK(myAlignment.maxloci,myAlignment.maxlocj);
		GAlign mydGAlign(strlen(seq1),strlen(seq2),params);
		mydGAlign.InitStrings(seq1,seq2,strlen(seq1),strlen(seq2));
		dinkelbach myDinkel(params, &mydGAlign);
		myDinkel.iteration(tempK); 
		answer =  mydGAlign.GetMeltingTempC(mydGAlign.maxloci,mydGAlign.maxlocj);
	}
	removeDollarSigns(seq1);
	removeDollarSigns(seq2);
	return answer;
}

int generateNonDegSet(char answer [][PRIMERS_MAX_SIZE], char * sequence) {
	char seed [PRIMERS_MAX_SIZE];
	int facDeg = 0;
	int lans = 1;

	for (int i=0; i<strlen(sequence);i++) {
		if (sequence[i] == 'N') {
			facDeg++;
			lans*=4;
		}
	}
	if(facDeg == 0) {
		strcpy (answer[0],sequence);
		return 1;
	} else if (lans >MAX_DEGENERACY) {
		return 0;
	}
	for (int h=0; h < lans;h++) {
		buildSeed(seed,h,facDeg);
		strcpy (answer[h],sequence);
		int seedIndex = 0;
		for (int i=0; i<strlen(sequence);i++) {
			if (sequence[i] == 'N') {
				answer [h][i] = seed[seedIndex];
				seedIndex++;
			}
		}
	}	
	return lans;

}

int buildCandidate (char * candidate, char * sequence, int finalPosition,int length) {
	strncpy(candidate,sequence+(finalPosition-length+1),length);
	candidate[length] = 0;
	int ndeg = strlen(degeneracy);
	for(int k=0;k<ndeg&&k<length;k++) {
		if(degeneracy[k] == '4' ) {
			candidate[length-k-1] = 'N';
		}
	}
	//Avoid candidates with N in the seed
	int nMask = strlen(maskTargets);
	for(int k=0;k<nMask && k<length;k++) {
		if (candidate[length-k-1] == 'N' && maskTargets [k] == '1') {
			return 0;
		}
	}
	return 1;
}

int calculateSize(char * mask) {
	int answer =0;
	for(int k=0; mask[k]!=0;k++) {
		if(mask[k]== '1' ) {
			answer++;
		}
	}
	return answer;
}

// Print methods
void   print_usage( char*  cmd )
{
  fprintf( stderr, "Usage: %s [OPTIONS]\n", cmd );
  fprintf( stderr, "OPTIONS:\n" );
  fprintf( stderr, "\t\t-tf FilePath\t\t\t: Targets file \n" );
  fprintf( stderr, "\t\t-nf FilePath\t\t\t: Non Targets file \n" );
  fprintf( stderr, "\t\t-pf FilePath\t\t\t: Primers file] \n" );
  fprintf( stderr, "\t\t-pof FilePath\t\t\t: Primers Output File \n" );
  fprintf( stderr, "\t\t-minPrimerLength N\t\t: Minimum Primer Length \n" );
  fprintf( stderr, "\t\t-maxPrimerLength N\t\t: Maximum Primer Length \n" );
  fprintf( stderr, "\t\t-minProdLength N\t\t: Minimum Product Length \n" );
  fprintf( stderr, "\t\t-maxProdLength N\t\t: Maximum Product Length \n" );
  fprintf( stderr, "\t\t-forwardPrimer S\t\t: Forward Primer \n" );
  fprintf( stderr, "\t\t-reversePrimer S\t\t: Reverse Primer \n" );
  fprintf( stderr, "\t\t-beginPosForward N\t\t: Begin of Range of Positions for Forward Primers \n" );
  fprintf( stderr, "\t\t-endPosForward N\t\t: End of Range of Positions for Forward Primers \n" );
  fprintf( stderr, "\t\t-beginPosReverse N\t\t: Begin of Range of Positions for Reverse Primers \n" );
  fprintf( stderr, "\t\t-endPosReverse N\t\t: End of Range of Positions for Reverse Primers \n" );
  fprintf( stderr, "\t\t-tmask S\t\t\t: Targets Mask \n" );
  fprintf( stderr, "\t\t-nmask S\t\t\t: Non Targets Mask \n" );
  fprintf( stderr, "\t\t-dmask S\t\t\t: Degeneracy Mask\n" );
  fprintf( stderr, "\t\t-numSourceSeq N\t\t\t: Number of Source Sequences \n" );
  fprintf( stderr, "\t\t-minCoverageTargets F\t\t: Minimum Targets Coverage \n" );
  fprintf( stderr, "\t\t-maxCoverageNonTargets F\t: Maximum Non Targets Coverage \n" );
  fprintf( stderr, "\t\t-maxSelfScore N\t\t\t: Maximum Self Complementarity Score \n" );
  fprintf( stderr, "\t\t-maxEndScore N\t\t\t: Maximum 3' End Complementarity Score \n" );
  fprintf( stderr, "\t\t-minGCContent F\t\t\t: Minimum GC Content \n" );
  fprintf( stderr, "\t\t-maxGCContent F\t\t\t: Maximum GC Content \n" );
  fprintf( stderr, "\t\t-gcClamp N\t\t\t: Minimum GC bases on the 3' end \n" );
  fprintf( stderr, "\t\t-maxPolyX N\t\t\t: Maximum mononucleotide repeat \n" );
  fprintf( stderr, "\t\t-primersConc F\t\t\t: Primer Concentration(M) \n" );
  fprintf( stderr, "\t\t-templateConc F\t\t\t: Template Concentration(M) \n" );
  fprintf( stderr, "\t\t-saltConc F\t\t\t: Salt Concentration(M) \n" );
  fprintf( stderr, "\t\t-saltCorrMethod N\t\t: Salt Correction Method. 1 for Santalucia's formula, 2 for Owczarzy's formula. \n" );
  fprintf( stderr, "\t\t-minTempTargets F\t\t: Targets Minimum Melting Temperature (C) \n" );
  fprintf( stderr, "\t\t-maxTempTargets F\t\t: Targets Maximum Melting Temperature (C) \n" );
  fprintf( stderr, "\t\t-maxTempNonTargets F\t\t: Non Targets Maximum Melting Temperature (C) \n" );
  fprintf( stderr, "\t\t-deltaTempNonTargets F\t\t: Non Targets Delta Melting Temperature(C) \n" );
  fprintf( stderr, "\t\t-maxPairTempDiff F\t\t: Max Pair Melting Temperature Difference(C) \n" );
  fprintf( stderr, "\t\t-primersLabel S\t\t\t: Primers Label Prefix \n" );
  fprintf( stderr, "\t\t-full_stats\t\t\t: Full Statistis on\n" );
  fprintf( stderr, "\t\t-h\t\t\t\t: Prints usage\n" );
  fflush(stdout);
}

void   printParameters(  ) {
  printf( "Parameters \n");
  printf( "Targets File: %s\n", targetsFile  );
  printf( "Non Targets File: %s\n", nonTargetsFile  );
  printf( "Primers File: %s\n", primersFile  );
  printf( "Primers Output File: %s\n", primersOutputFile  );
  printf( "Min Primer Length: %d\n", minLength  );
  printf( "Max Primer Length: %d\n", maxLength  );
  printf( "Min Product Length: %d\n", productMinLength  );
  printf( "Max Product Length: %d\n", productMaxLength  );
  printf( "Forward Primer: %s\n", forwardPrimer  );
  printf( "Reverse Primer: %s\n", reversePrimer  );
  printf( "Range of positions for forward primers: %d - %d\n", beginPosForward, endPosForward  );
  printf( "Range of positions for reverse primers: %d - %d\n", beginPosReverse, endPosReverse  );
  printf( "Targets Mask: %s\n", maskTargets  );
  printf( "Non Targets Mask: %s\n", maskNonTargets  );
  printf( "Degeneracy Mask: %s\n", degeneracy  );
  printf( "Number of Source Sequences: %d\n", nseqForCandidates );
  printf( "Min Coverage Targets: %.2f\n", minCoverageTargets );
  printf( "Max Coverage Non Targets: %.2f\n", maxCoverageNonTargets );
  printf( "Max Self Score: %.2f\n", maxSelfScore  );
  printf( "Max End Score: %.2f\n", maxEndScore  );
  printf( "Min GC Content: %.2f\n", minGCContent  );
  printf( "Max GC Content: %.2f\n", maxGCContent  );
  printf( "GC Clamp: %d\n", gcClamp  );
  printf( "Max Poly-X: %d\n", maxPolyX  );
  printf( "Primer Concentration: %.9f\n", concPrimers  );
  printf( "Template Concentration: %.9f\n", concSequences  );
  printf( "Salt Concentration: %.5f\n", salt  );
  printf( "Salt Correction Method: %d\n", salMethod  );
  printf( "Targets Min Melting Temp : %.2f\n", minTempTargets  );
  printf( "Targets Max Melting Temp: %.2f\n", maxTempTargets  );
  printf( "Non Targets Max Melting Temp: %.2f\n", maxTempNonTargets  );
  printf( "Non Targets Delta Temp: %.2f\n", deltaTempNonTargets  );
  printf( "Max Pair Temperature Difference: %.2f\n", maxPairTempDiff  );
  printf( "Primers Label: %s\n", primersLabel  );
  if (fullStats) {
	printf( "Full stats\n" );
  } else {
	printf( "Dependent stats\n" );
  }
  printf( "Sequence for ranking: %s\n", seqRank  );
  fflush(stdout);
  
}

void printTargets() {
	printf("Targets: %d\n",nTargets);
	for (int i=0;i<nTargets;i++) {
		printf("%s\n",targets[i]);
	}
}
void printSeedsInfo () {
	int total = 0;
	printf("Seeds: %d\n",nSeedTargets);
	for (int i=0;i<nSeedTargets;i++) {
		OccuranceList list = (seedsTableTargets [i]).occurancesList;
		if( list->size >0 ) {
			printf("Seed: %s, Occurances: %d\n",(seedsTableTargets [i]).seed,list->size);
			for (Occurance ptr = list->head;ptr!=NULL;ptr = ptr->next) {
				printf("%d-%d ",ptr->sequenceId+1,ptr->sequencePosition);
			}
			total += list->size;
			printf("\n");	
		}
	}
	printf("Total occurances: %d\n",total);	
}


void printNonTargetsSeedsInfo () {
	int total = 0;
	printf("Seeds: %d\n",nSeedNonTargets);
	for (int i=0;i<nSeedNonTargets;i++) {
		OccuranceList list = (seedsTableNonTargets [i]).occurancesList;
		if( list->size >0 ) {
			printf("Seed: %s, Occurances: %d\n",(seedsTableNonTargets [i]).seed,list->size);
			for (Occurance ptr = list->head;i<1&&ptr!=NULL;ptr = ptr->next) {
				printf("%d-%d ",ptr->sequenceId+1,ptr->sequencePosition);
			}
			total += list->size;
			printf("\n");	
		}
	}
	printf("Total occurances: %d\n",total);	
}

float calculatePercentageTargetsCovered(Primer primer) {
	float percentage = 0;
	for(int i=0;i<primer->nPositions;i++) {
		if( primer->statusTargets[i]==1) {
			percentage++;
		}
	}
	if(primer->nPositions > 0 ) {
		percentage = percentage*100 /(float)primer->nPositions;
	}
	return percentage;
}

float calculatePercentageNonTargetsCovered(Primer primer) {
	float percentage = 0;
	for(int i=0;i<primer->nNonTargets;i++) {
		if( primer->statusNonTargets[i]==1) {
			percentage++;
		}
	}
	if(primer->nNonTargets > 0 ) {
		percentage = percentage*100 /(float)primer->nNonTargets;
	}
	return percentage;
}

void printPrimer (Primer p) {
	char seqName [10];
	if(p->seqOrig>=0) {
		sprintf(seqName,"S%d",p->seqOrig+1);
	} else {
		sprintf(seqName,"G");
	}
	float perTC = calculatePercentageTargetsCovered(p);
	float perNTC = calculatePercentageNonTargetsCovered(p);
	printf("%s_%s_%s_%d ", p->forward?"Forward":"Reverse",primersLabel,seqName,p->id);
	printf("%s %d %.2f %.2f %.2f %.2f %.2f%% %.2f%%",p->seq,p->positions[seqPositions],p->minTempTargets,p->avgTempTargets,p->maxTempTargets,p->maxTempNonTargets, perTC, perNTC);
	if(strlen(seqRank) >0) {
		printf("%.2f %.2f%%",p->tempRanking, 100*p->ranking);
	}
}

void printPrimers() {
	for (int i=0;i<nPrimers;i++) {
		printPrimer(primers[i]);
		printf("\n");
	}
	printf("Primers: %d\n", nPrimers);
}

float calculatePercentageTargetsCoveredByPair(Pair p) {
	float percentage = 0;
	Primer p1 = p->forwardPrimer;
	Primer p2 = p->reversePrimer;
	for(int i=0;i<p1->nPositions;i++) {
		if( p1->statusTargets[i]==1 && p2->statusTargets[i]==1) {
			percentage++;
		}
	}
	if(p1->nPositions > 0 ) {
		percentage = percentage*100 /(float)p1->nPositions;
	}
	return percentage;
}

float calculatePercentageNonTargetsCoveredByPair(Pair p) {
	float percentage = 0;
	Primer p1 = p->forwardPrimer;
	Primer p2 = p->reversePrimer;
	for(int i=0;i<p1->nNonTargets;i++) {
		if( p1->statusNonTargets[i]==1 && p2->statusNonTargets[i]==1) {
			percentage++;
		}
	}
	if(p1->nNonTargets > 0 ) {
		percentage = percentage*100 /(float)p1->nNonTargets;
	}
	return percentage;
}


void printPairs() {
	for (int i=0;i<nPairs;i++) {
		printPrimer(pairs[i] -> forwardPrimer);
		printf(" - ");
		printPrimer(pairs[i] -> reversePrimer);
		float perTC = calculatePercentageTargetsCoveredByPair(pairs[i]);
		float perNTC = calculatePercentageNonTargetsCoveredByPair(pairs[i]);
		printf(" - %.2f%% %.2f%%", perTC, perNTC);
		printf("\n");
	}
	printf("Pairs: %d\n", nPairs);
	fflush(stdout);
}

void printStatistics (Statistics stats) {
	double total = stats -> counts [STAT_TOTAL_CANDIDATES];
	double percentages [STAT_COUNTS];
	for (int i = 0 ; i<STAT_COUNTS;i++ ) {
		percentages[i] = 100*((double) stats ->counts [i])/total;
	}
	printf("Statistics\n");
	printf("Total candidates: %d\n",stats -> counts [STAT_TOTAL_CANDIDATES]);
	printf("Pass degeneracy %d (%.2f \%)\n", stats -> counts [STAT_TEST_DEGENERACY], percentages[STAT_TEST_DEGENERACY]);
	printf("Pass GC Content %d (%.2f \%)\n", stats -> counts [STAT_TEST_GCCONTENT], percentages[STAT_TEST_GCCONTENT]);
	printf("Pass GC Clamp %d (%.2f \%)\n", stats -> counts [STAT_TEST_GCCLAMP], percentages[STAT_TEST_GCCLAMP]);
	printf("Pass Max Poly-X %d (%.2f \%)\n", stats -> counts [STAT_TEST_MAX_POLY_X], percentages[STAT_TEST_MAX_POLY_X]);
	printf("Pass self complementarity %d (%.2f \%)\n", stats -> counts [STAT_TEST_SELFCOMP], percentages[STAT_TEST_SELFCOMP]);
	printf("Pass self melting %d (%.2f \%)\n", stats -> counts [STAT_TEST_TMSELF], percentages[STAT_TEST_TMSELF]);
	printf("\tFailed self melting low %d (%.2f \%)\n", stats -> counts [STAT_FAILED_SELFCOMP_LOW], percentages[STAT_FAILED_SELFCOMP_LOW]);
	printf("\tFailed self melting high %d (%.2f \%)\n", stats -> counts [STAT_FAILED_SELFCOMP_HIGH], percentages[STAT_FAILED_SELFCOMP_HIGH]);
	printf("Pass targets %d (%.2f \%)\n", stats -> counts [STAT_TEST_TARGETS], percentages[STAT_TEST_TARGETS]);
	printf("Pass non targets %d (%.2f \%)\n", stats -> counts [STAT_TEST_NONTARGETS], percentages[STAT_TEST_NONTARGETS]);
	printf("Pass all tests %d (%.2f \%)\n\n", stats -> counts [STAT_ALLTESTS], percentages[STAT_ALLTESTS]);
	/*printf("Targets test fail details\n");
	for(int i = 0; i < nTargets; i++) {
		printf("Target: %d throw %d\n",i,stats -> failTargetsCounts [i]);
	}*/
	fflush(stdout);
}

//List management methods

void addOccurance (OccuranceList occurancesList, int sequenceId, int sequencePosition ) {
	Occurance oc = (Occurance)malloc (sizeof(TOccurance));
	oc -> sequenceId = sequenceId; 
	oc -> sequencePosition = sequencePosition;
	oc -> next = NULL;
	if(occurancesList -> head == NULL) {
		occurancesList -> head = oc;
	}
	if(occurancesList -> tail != NULL) {
		occurancesList -> tail -> next = oc;
	}
	occurancesList -> tail = oc;
	(occurancesList -> size)++;
}

void deleteList (Occurance occurancesList) {
	if(occurancesList==NULL) {
		return;
	}
	deleteList(occurancesList->next);
	free (occurancesList);
}

//Solution methods

void initialize() {
	for (int i=0;i<SEEDSTABLE_TAM;i++) {
		(seedsTableTargets [i]).occurancesList = (OccuranceList)malloc (sizeof(TOccuranceList));
		(seedsTableTargets [i]).occurancesList -> head = (seedsTableTargets [i]).occurancesList -> tail = NULL;
		(seedsTableTargets [i]).occurancesList -> size = 0;
	}
	for (int i=0;i<SEEDSTABLE_TAM;i++) {
		(seedsTableNonTargets [i]).occurancesList = (OccuranceList)malloc (sizeof(TOccuranceList));
		(seedsTableNonTargets [i]).occurancesList -> head = (seedsTableNonTargets [i]).occurancesList -> tail = NULL;
		(seedsTableNonTargets [i]).occurancesList -> size = 0;
	}
}

void depureSequence(char * sequence) {
	int n = strlen(sequence);
	int j=0;
	for (int i=0; i<n;i++) {
		if(sequence[i]!=' '&& !iscntrl(sequence[i])) {
			sequence[j] = toupper(sequence[i]);
			if (sequence[j] != 'A' && sequence[j] != 'C' && sequence[j] != 'G' && sequence[j] != 'T' && sequence[j] != 'N') {
				sequence[j] = 'N';
			}
			j++;
		}
	}
	sequence[j] = 0;
}

void loadFile (char * fileName,int isTarget, int sequences) {
	char line [SEQUENCES_MAX_SIZE];
	char nextSequence [SEQUENCES_MAX_SIZE];
	nextSequence[0] = 0;
	FILE * file = fopen (fileName, "rt");
	if (file == NULL) {
		printf("File %d could not be opened\n",fileName);
		return;
	}
	fgets(line,SEQUENCES_MAX_SIZE,file);
	int j;
	for(j=0;j<sequences&&!feof(file);) {
		if(line[0]=='>') {
			depureSequence(nextSequence);
			if(strlen(nextSequence) > 0) {
				if(isTarget) {
					targets[nTargets] = (char *)malloc ((strlen(nextSequence)+1)*sizeof(char));
					strcpy(targets[nTargets],nextSequence);
					nTargets++;
					if(nTargets == TARGETS_MAX_SIZE) {
						printf("Maximum supported number of targets reached\n");
						break;
					}
				} else {
					nontargets[nNonTargets] = (char *)malloc ((strlen(nextSequence)+1)*sizeof(char));
					strcpy(nontargets[nNonTargets],nextSequence);
					nNonTargets++;
					if(nNonTargets == NON_TARGETS_MAX_SIZE) {
						printf("Maximum supported number of non targets reached\n");
						break;
					}
				}
				j++;
				nextSequence[0] = 0;
			}
			fgets(line,SEQUENCES_MAX_SIZE,file);
		}
		strcat(nextSequence,line);
		fgets(line,SEQUENCES_MAX_SIZE,file);
	}
	fclose(file);
	depureSequence(nextSequence);
	if(j<sequences && strlen(nextSequence) > 0) {
		if(isTarget && nTargets < TARGETS_MAX_SIZE) {
			targets[nTargets] = (char *)malloc ((strlen(nextSequence)+1)*sizeof(char));
			strcpy(targets[nTargets],nextSequence);
			nTargets++;
		} else if (!isTarget && nNonTargets < NON_TARGETS_MAX_SIZE){
			nontargets[nNonTargets] = (char *)malloc ((strlen(nextSequence)+1)*sizeof(char));
			strcpy(nontargets[nNonTargets],nextSequence);
			nNonTargets++;
		}
	}	 
}

void loadFiles (int sequences) {
	if(strlen(targetsFile)>0) {
		loadFile(targetsFile,1,sequences);
	}
	if(strlen(nonTargetsFile)>0) {
		loadFile(nonTargetsFile,0,sequences);
	}
}

void addSequenceInfo(char * sequence, char * mask, TSeed * seedsTable, int sequenceId ) {
	char seed [SEEDS_MAX_SIZE+1];
	int nmask = strlen(mask);
	if (nmask == 0) {
		return;
	}
	int nseq = strlen(sequence);
	for(int i=nmask-1;i<nseq;i++) {
		int seedIndex = 0;
		for(int j=0;j<nmask;j++) {
			if (mask[j] == '1') {
				if(sequence[i-j]=='N') {
					seedIndex = -1;
					break;
				}
				seed[seedIndex] = sequence[i-j];	
				seedIndex++;
			}
		}
		if(seedIndex >=0) {
			seed[seedIndex]=0;
			int hash = getHash(seed);
			if ((seedsTable[hash]).occurancesList -> size == 0) {
				strcpy((seedsTable[hash]).seed,seed);
			}
			addOccurance ((seedsTable[hash]).occurancesList, sequenceId, i );
		}
		
	}
}

void depureSeedsTargets(int nseq) {
	
	for(int i=0; i<nSeedTargets;i++) {
		int deleteSeed = 0;
		OccuranceList list = (seedsTableTargets [i]).occurancesList;
		if (list->size < nseq) {
			deleteSeed = 1;
		} else {
			int previousSequence = -1;
			for (Occurance ptr = list->head;ptr!=NULL;ptr = ptr->next) {
				if(ptr->sequenceId -previousSequence > 1) {
					deleteSeed = 1;
					break;
				}
				previousSequence = ptr->sequenceId;
			}
			if(previousSequence < nseq-1) {
				deleteSeed = 1;
			}
		}
		if(deleteSeed) {
			deleteList(list->head);
			list -> size = 0;
			list -> head = list -> tail = NULL; 
		}
	}
}
void buildSeedsTable (char * sequences [], int nseq, char * mask, int areTargets) {
	int size = calculateSize(mask);
	int size2 = (int)pow((double)4, (double)size);
	TSeed * seedsTable;
	if(areTargets) {
		seedsTable = seedsTableTargets;
		nSeedTargets = size2;
	} else {
		seedsTable = seedsTableNonTargets;
		nSeedNonTargets = size2;
	}

	for(int i=0; i<nseq;i++) {
		addSequenceInfo(sequences[i],mask,seedsTable,i);
	}
}

void buildTables () {
	buildSeedsTable(targets, nTargets, maskTargets, 1);
	buildSeedsTable(nontargets, nNonTargets, maskNonTargets, 0);
	depureSeedsTargets(nTargets);	
}



void set_dpal_args(dpal_args * a) {
    unsigned int i, j;        

    memset(a, 0, sizeof(*a));
    for (i = 0; i <= UCHAR_MAX; i++)
        for (j = 0; j <= UCHAR_MAX; j++)
          if (('A' == i || 'C' == i || 'G' == i || 'T' == i || 'N' == i)
              && ('A' == j || 'C' == j || 'G' == j || 'T' == j
                  || 'N' == j)) {
            if (i == 'N' || j == 'N')
              a->ssm[i][j] = -25;      
            else if (i == j)
              a->ssm[i][j] = 100;
            else         
              a->ssm[i][j] = -100; 
          } else
            a->ssm[i][j] = INT_MIN;

    a->gap                = -200;
    a->gapl               = -200;
    a->flag               = DPAL_LOCAL;
    a->max_gap            = 1;
    a->fail_stop          = 1;
    a->check_chars        = 0;
    a->debug              = 0;
    a->score_only         = 1;
    a->force_generic      = 0;
    a->force_long_generic = 0;
    a->force_long_maxgap1 = 0;
}

int calculatePosition (int seqLength, int index, int forward) {
	if(forward) {
		return index;
	} else {
		return seqLength-index-1;
	}
}

int testDegeneracy ( char * candidate) {
	int lseq = strlen(candidate);
	int ldeg = strlen(degeneracy);
	//Check degeneracy
	for(int k=0;k<lseq;k++) {
		if (candidate[lseq-k-1] == 'N' && (k >= ldeg || degeneracy[k] != '4') ) {
			return 0;
		}
	}
	return 1;
}

int testGCContent ( char * candidate) {
	int lseq = strlen(candidate);
	float count = countGCContent(candidate);
	float mean = 100*count/(float)lseq;
	if(mean>=minGCContent && mean <= maxGCContent) {
		return 1;
	}
	return 0;
}

int testGCClamp ( char * candidate) {
	int lseq = strlen(candidate);
	for(int k=0;k<lseq && k<gcClamp;k++) {
		if (candidate[lseq-k-1] != 'C' && candidate[lseq-k-1] != 'G' ) {
			return 0;
		}
	}
	return 1;
}

int testMaxPolyX ( char * candidate) 
{
	int lseq = strlen(candidate);
	int repeats=0;
	char lastBase = ' ';
	for(int k=0;k<lseq ;k++) 
	{
		if(candidate[k] == lastBase) 
		{
			repeats++;
			if (repeats > maxPolyX) 
			{
				return 0;
			}
		}
		else 
		{
			repeats = 1;
		}
		lastBase = candidate[k];
	}
	return 1;
}

int testComplementarity (char nonDegCandidate [][PRIMERS_MAX_SIZE], int lenSetCand, char sequence[]) {
	int lseq = strlen(nonDegCandidate[0]);
	int lseq2 = strlen(sequence);
	char reverseComplement [PRIMERS_MAX_SIZE];
	char nonDegSequence [MAX_DEGENERACY][PRIMERS_MAX_SIZE];
	int lenSetSeq = generateNonDegSet (nonDegSequence,sequence);
	if(lenSetSeq == 0) {
		return 0;
	}
	dpal_results r;
	for(int i = 0; i < lenSetSeq; i++ ) {
		calculateReverseComplement( reverseComplement, nonDegSequence[i] );
        	for(int j = 0; j < lenSetCand; j++ ) {
			// Test of self complementarity adapted from primer3 
			_dpal_long_nopath_maxgap1_local(nonDegCandidate[j], reverseComplement, lseq, lseq2, &local_args, &r);
			if( r.score > maxSelfScore ) {
				return 0;
			}
			_dpal_long_nopath_maxgap1_global_end(nonDegCandidate[j], reverseComplement, lseq, lseq2, &end_args, &r);
			if( r.score > maxEndScore ) {
				return 0;
			}

		}
	}
	return 1;
}
int testMeltingTemperature (char nonDegCandidate [][PRIMERS_MAX_SIZE], int lenSetCand, Statistics stats) {
	char complement [PRIMERS_MAX_SIZE];
        for(int j = 0; j < lenSetCand; j++ ) {
		char * seq = nonDegCandidate[j];
		calculateComplement(complement,seq);
		// check primer Tm is between minTempTargets and maxTempTargets
		float melting = calculateMeltingTemperature(seq, complement,paramsPrimerSeq);
		if (melting<minTempTargets ){
			stats->counts [STAT_FAILED_SELFCOMP_LOW] += 1;
			return 0;
		}
		if (melting>maxTempTargets) {
			stats->counts [STAT_FAILED_SELFCOMP_HIGH] += 1;
			return 0;
		}
	}
	return 1;
}

void updateTemperatureData(Primer primer,float temp, int targets) {
	if(targets) {
		if(primer -> minTempTargets > temp) {
			primer->minTempTargets = temp;
		}
		if(primer -> maxTempTargets < temp) {
			primer->maxTempTargets = temp;
		}
	} else {
		if(primer -> minTempNonTargets > temp) {
			primer->minTempNonTargets = temp;
		}
		if(primer -> maxTempNonTargets < temp) {
			primer->maxTempNonTargets = temp;
		}
	}
} 
int checkMeltingTemperatureTargets(char nonDegCandidate [][PRIMERS_MAX_SIZE], int lenSetCand, char * sequence, float * outTempLocal) {
	char nonDegSequence [MAX_DEGENERACY][PRIMERS_MAX_SIZE];
	int lenSetSeq = generateNonDegSet (nonDegSequence,sequence);
	
	if(lenSetSeq == 0) {
		*outTempLocal = -274;
		return 0;
	}
	int answer = 1;
	*outTempLocal = maxTempTargets+1;
	for(int i = 0; i < lenSetSeq; i++ ) {
		int seqCovered = 0;
		float maxTempNonDegSeq = -274;
		for(int j = 0; j < lenSetCand; j++ ) {
			float melting = calculateMeltingTemperature(nonDegCandidate [j],nonDegSequence[i],paramsPrimerSeq);
			if(maxTempNonDegSeq < melting) {
				maxTempNonDegSeq = melting;
			}
			if (melting>=minTempTargets && melting<=maxTempTargets) {
				seqCovered = 1;
			}
		}
		if (*outTempLocal > maxTempNonDegSeq ) {
			*outTempLocal = maxTempNonDegSeq;
		}
		if(!seqCovered) {
			answer = 0;
		}
	}
	return answer;	
}

void updateAverageTemperature(Primer primer) {
	float average = 0;
	for(int i=0;i<primer->nPositions;i++) {
		float nextTemp = primer->tempsTargets[i];
		average+=nextTemp;
	}
	if(primer->nPositions > 0 ) {
		average = average /(float)primer->nPositions;
		primer->avgTempTargets = average;
	}
}
void calculateSequenceForMelting (char answer [], char sequence [], int finalPosition, int length, int extraBases) {
	int seqLen = strlen(sequence);
	answer[0] = 0;
	//Calculating initial position without extra bases
	int initialPosition = finalPosition-length+1;
	//Add extra bases
	initialPosition-=extraBases;
	finalPosition+=extraBases;
	//Fix if limits are exceeded 
	if(initialPosition<0) {
		initialPosition = 0;
	}
	if(finalPosition > seqLen-1) {
		finalPosition = seqLen-1;
	}
	if(finalPosition < initialPosition ) {
		return;
	}
	//Recalculate length
	length = finalPosition-initialPosition+1;
	strncpy(answer,sequence+initialPosition,length);
	answer[length]=0;
	//Calculate complement
	char tmp [2*PRIMERS_MAX_SIZE];
	calculateComplement(tmp,answer);
	strcpy(answer,tmp);
}   
int testCandidateAgainstTargets(char nonDegCandidate [][PRIMERS_MAX_SIZE], int lenSetCand, Primer primer, Statistics stats, int verbose ) {
	char nextToCompare [2*PRIMERS_MAX_SIZE];
	char seed [SEEDS_MAX_SIZE+1];
	float maxTempSeq = minTempTargets-1;
	int maxSeqNonCoveredAllowed = (int)((float)(100-minCoverageTargets)*(float)nTargets/(float)100); 
	primer -> minTempTargets = maxTempTargets+1;
	primer -> maxTempTargets = minTempTargets-1;
	//All non deg candidates must share the seed
	char * candidate = nonDegCandidate [0];
	buildSeedFromSequence(candidate, maskTargets, seed);
	int index = getHash(seed);
	int length = strlen(candidate);
	int sequenceId=-1;
	int lastCoveredSeqId=-1;
	int targetsNonCovered=0;
	OccuranceList list = (seedsTableTargets[index]).occurancesList;
	if(list == NULL || list->size == 0) {
		return 0;
	}
	for (Occurance ptr = list->head;ptr!=NULL;ptr = ptr->next) {
		if(ptr->sequenceId > lastCoveredSeqId + 1) {
			targetsNonCovered += (ptr->sequenceId-lastCoveredSeqId-1);
			if (verbose) printf("Candidate %s does not cover target %d\n",primer->seq,ptr->sequenceId);
			(stats->failTargetsCounts[sequenceId+1])++;
			lastCoveredSeqId = ptr->sequenceId - 1; 
			if(targetsNonCovered > maxSeqNonCoveredAllowed) {
				return 0;
			}
		}
		char * sequence = targets[ptr->sequenceId];
		int seqLen = strlen(sequence);
		int position = ptr->sequencePosition;
		//Two extra bases added on each side to avoid ignoring dangling ends effect
		calculateSequenceForMelting(nextToCompare,sequence,position,length,2);
		float melting;
		int pass = checkMeltingTemperatureTargets(nonDegCandidate,lenSetCand,nextToCompare,&melting);
		int hybPos = calculatePosition (seqLen, position, primer -> forward);
		
		//if(verbose) printf("Candidate %s reports %f for target %d in position %d\n",primer->seq,melting,ptr->sequenceId+1,hybPos);
		if(ptr->sequenceId == sequenceId ) {
			if(melting>maxTempSeq) {
				maxTempSeq=melting;
				primer->positions[sequenceId] = hybPos;
				primer->tempsTargets[sequenceId] = melting;
			}
		} else {
			if(sequenceId >=0) {
				//Compare temperature of the previous sequence against global maximum and minimum
				updateTemperatureData(primer,maxTempSeq,1);
				if(verbose) printf("Maximum temperature for candidate %s in target %d is %.2f in position %d\n",primer->seq,sequenceId+1,primer->tempsTargets[sequenceId],primer->positions[sequenceId]);
			} 
			//Initialize temperature for the new sequence
			maxTempSeq = melting;
			sequenceId = ptr->sequenceId;
			primer->positions[sequenceId] = hybPos;
			primer->tempsTargets[sequenceId] = melting;
		}
		if(pass) {
			lastCoveredSeqId = ptr->sequenceId;
			primer->statusTargets[ptr->sequenceId] = 1;
		}
	}
	if(lastCoveredSeqId < nTargets-1) {
		targetsNonCovered += (nTargets-lastCoveredSeqId-1);
		if(verbose) printf("Candidate %s does not cover target %d\n",primer->seq,nTargets);
		(stats->failTargetsCounts[sequenceId+1])++;
		if(targetsNonCovered > maxSeqNonCoveredAllowed) {
			return 0;
		}	
	}
	//Compare temperature of the last sequence against global maximum and minimum
	updateTemperatureData(primer,maxTempSeq,1);
	if(verbose) printf("Maximum temperature for candidate %s in target %d is %.2f in position %d\n",primer->seq,sequenceId+1,primer->tempsTargets[sequenceId],primer->positions[sequenceId]);
	updateAverageTemperature(primer);
	return 1;
}

float calculateMinTempCovered(Primer primer) {
	float minTempCovered = 10000;
	for(int i=0;i<primer->nPositions;i++) {
		if(primer->statusTargets[i]) {
			float melting = primer->tempsTargets[i];
			if (melting<minTempCovered) {
				minTempCovered = melting;
			}
		}
	}
	return minTempCovered;
} 
int checkMeltingTemperatureNonTargets(char nonDegCandidate [][PRIMERS_MAX_SIZE], int lenSetCand, char * sequence, Primer primer,int nonTargetId) {
	char nonDegSequence [MAX_DEGENERACY][PRIMERS_MAX_SIZE];
	int lenSetSeq = generateNonDegSet (nonDegSequence,sequence);
	if(lenSetSeq == 0) {
		return 0;
	}
	float minForDelta = calculateMinTempCovered(primer);
	for(int i = 0; i < lenSetSeq; i++ ) {
		for(int j = 0; j < lenSetCand; j++ ) {
			float melting = calculateMeltingTemperature(nonDegCandidate [j],nonDegSequence[i],paramsPrimerSeq);
			if (primer->tempsNonTargets[nonTargetId]< melting) {
				primer->tempsNonTargets[nonTargetId]=melting;
			}
			updateTemperatureData(primer,melting,0);
			if(melting > maxTempNonTargets) {
				return 0;
			}
			if (melting > minForDelta - deltaTempNonTargets) {
				return 0;
			}
		}
	}
	return 1;	
}


int testCandidateAgainstNonTargets( char nonDegCandidate [][PRIMERS_MAX_SIZE], int lenSetCand, Primer primer) {
	if (nNonTargets == 0) {
		return 1;
	}
	char nextToCompare [2*PRIMERS_MAX_SIZE];
	char seed [SEEDS_MAX_SIZE+1];
	int maxSeqCoveredAllowed = (int)((float)(maxCoverageNonTargets)*(float)nNonTargets/(float)100); 
	primer -> minTempNonTargets = 10000;
	primer -> maxTempNonTargets = -274;

	//All non deg candidates must share the seed
	char * candidate = nonDegCandidate [0];
	buildSeedFromSequence(candidate, maskNonTargets, seed);
	int index = getHash(seed);
	int length = strlen(candidate);
	int nonTargetsCovered=0;
	OccuranceList list = (seedsTableNonTargets[index]).occurancesList;
	if(list == NULL || list->size == 0) {
		return 1;
	}
	int lastSequenceId = -1;
	for (Occurance ptr = list->head;ptr!=NULL;ptr = ptr->next) {
		if(lastSequenceId < ptr->sequenceId) {
			primer->tempsNonTargets[ptr->sequenceId]=-274;
			lastSequenceId = ptr->sequenceId;
		}
		char * sequence = nontargets[ptr->sequenceId];
		int position = ptr->sequencePosition;
		//Two extra bases added on each side to avoid ignoring dangling ends effect
		calculateSequenceForMelting(nextToCompare,sequence,position,length,2);
		if(!checkMeltingTemperatureNonTargets(nonDegCandidate,lenSetCand,nextToCompare,primer,ptr->sequenceId)) {
			nonTargetsCovered++;
			primer->statusNonTargets[ptr->sequenceId]=1;
			if(nonTargetsCovered > maxSeqCoveredAllowed) {
				return 0;
			}
		}
	}	
	return 1;
}

int testCandidateAgainstNonTargetsFull(char nonDegCandidate [][PRIMERS_MAX_SIZE], int lenSetCand, Primer primer) {
	primer -> minTempNonTargets = 10000;
	primer -> maxTempNonTargets = -274;
	if (nNonTargets == 0) {
		return 1;
	}
	int nonTargetsCovered=0;
	float minForDelta = calculateMinTempCovered(primer);
	char complement [SEQUENCES_MAX_SIZE];
	int maxSeqCoveredAllowed = (int)((float)(maxCoverageNonTargets)*(float)nNonTargets/(float)100); 

	for(int	i=0;i<nNonTargets;i++) {
		primer->tempsNonTargets[i]=-274;
		calculateComplement(complement,nontargets[i]);
		float melting;
		for(int j = 0; j < lenSetCand; j++ ) {
			melting = calculateMeltingTemperature(nonDegCandidate [j], complement,paramsPrimerSeq);
			if (primer->tempsNonTargets[i] < melting) {
				primer->tempsNonTargets[i]=melting;
			}
		}
		melting =primer->tempsNonTargets[i];
		updateTemperatureData(primer,melting,0);
		if (melting > maxTempNonTargets || melting>minForDelta-deltaTempNonTargets) {
			nonTargetsCovered++;
			primer->statusNonTargets[i]=1;
			if(nonTargetsCovered > maxSeqCoveredAllowed) {
				return 0;
			}
		}
	}	
	return 1;
}

int testLoadedCandidateAgainstTargets(Primer primer) {
	int maxSeqNonCoveredAllowed = (int)((float)(100-minCoverageTargets)*(float)primer->nPositions/(float)100); 
	int targetsNonCovered=0;
	for(int i=0;i<primer->nPositions;i++) {
		float melting = primer->tempsTargets[i];
		if (melting<minTempTargets || melting>maxTempTargets) {
			targetsNonCovered++;
			if(targetsNonCovered > maxSeqNonCoveredAllowed) {
				return 0;
			}
		} else {
			primer->statusTargets[i]=1;
		}
	}
	return 1;
}

int testLoadedCandidateAgainstNonTargets(Primer primer) {
	float minForDelta = calculateMinTempCovered(primer);
	int maxNonTargetsCoveredAllowed = (int)((float)(maxCoverageNonTargets)*(float)primer->nNonTargets/(float)100); 
	int nonTargetsCovered=0;
	for(int i=0;i<primer->nNonTargets;i++) {
		float melting = primer->tempsNonTargets[i];
		if (melting > maxTempNonTargets || melting>minForDelta-deltaTempNonTargets) {
			nonTargetsCovered++;
			primer->statusNonTargets[i]=1;
			if(nonTargetsCovered > maxNonTargetsCoveredAllowed) {
				return 0;
			}
		}
	}
	return 1;
}

void rankPrimer(char nonDegCandidate [][PRIMERS_MAX_SIZE], int lenSetCand, Primer primer) {
	float maxTemp = -274;
	char complement [SEQUENCES_MAX_SIZE];
	calculateComplement(complement,seqRank);
	for(int j = 0; j < lenSetCand; j++ ) {
		float melting = calculateMeltingTemperature(nonDegCandidate [j],seqRank,paramsPrimerSeq);
		if(maxTemp < melting) {
			maxTemp = melting;
		}
	}
	primer->tempRanking = maxTemp;
	float ranking = 0;
	for(int i=0;i<primer->nPositions;i++) {
		float nextTemp = primer->tempsTargets[i];
		if(nextTemp <= maxTemp) {
			ranking++;
		}
	}
	ranking = ranking/(float)primer->nPositions;
	primer->ranking = ranking;
}

int testCandidate ( Primer primer, Statistics stats, int verbose) {
	char nonDegCandidate [MAX_DEGENERACY][PRIMERS_MAX_SIZE];
	int lenSetCand = generateNonDegSet (nonDegCandidate,primer->seq);
	stats->counts [STAT_TOTAL_CANDIDATES] +=1;
	int passTests = 0;
	int onePass = testDegeneracy (primer->seq);
	stats->counts [STAT_TEST_DEGENERACY] += onePass;
	passTests += onePass;
	if(verbose && !onePass) printf("Candidate %s failed degeneracy test\n",primer->seq);
	if (fullStats || passTests == 1) {
		onePass = testGCContent(primer->seq);
		stats->counts [STAT_TEST_GCCONTENT] += onePass;	
		passTests += onePass;
		if(verbose && !onePass) printf("Candidate %s failed GC Content test\n",primer->seq);
	}
	if (fullStats || passTests == 2) {
		onePass = testGCClamp(primer->seq);
		stats->counts [STAT_TEST_GCCLAMP] += onePass;	
		passTests += onePass;
		if(verbose && !onePass) printf("Candidate %s failed GC Clamp test\n",primer->seq);
	}
	if (fullStats || passTests == 3) {
		onePass = testMaxPolyX(primer->seq);
		stats->counts [STAT_TEST_MAX_POLY_X] += onePass;	
		passTests += onePass;
		if(verbose && !onePass) printf("Candidate %s failed Max Poly-X test\n",primer->seq);
	}
	if (fullStats || passTests == 4) {
		onePass = testComplementarity(nonDegCandidate,lenSetCand,primer->seq);
		stats->counts [STAT_TEST_SELFCOMP] += onePass;	
		passTests += onePass;
		if(verbose && !onePass) printf("Candidate %s failed self complementarity test\n",primer->seq);
	}
	
	if (fullStats || passTests == 5) {
		onePass = testMeltingTemperature(nonDegCandidate,lenSetCand,stats);
		stats->counts [STAT_TEST_TMSELF] += onePass;	
		passTests += onePass;
		if(verbose && !onePass) printf("Candidate %s failed melting temperature test\n",primer->seq);
	}
	if (fullStats || passTests == 6) {
		if(primer->origin == 'L') {
			onePass = testLoadedCandidateAgainstTargets(primer);
		} else {
			onePass = testCandidateAgainstTargets(nonDegCandidate,lenSetCand,primer,stats,verbose);
		} 
		stats->counts [STAT_TEST_TARGETS] += onePass;	
		passTests += onePass;
		if(verbose && !onePass) printf("Candidate %s failed targets test\n",primer->seq);
	}
	if (fullStats || passTests == 7) {
		if(primer->origin == 'L') {
			onePass = testLoadedCandidateAgainstNonTargets(primer);
		} else if(strlen(maskNonTargets)>0) {
			onePass = testCandidateAgainstNonTargets(nonDegCandidate,lenSetCand,primer);
		} else {
			onePass = testCandidateAgainstNonTargetsFull(nonDegCandidate,lenSetCand,primer);
		}
		stats->counts [STAT_TEST_NONTARGETS] += onePass;	
		passTests += onePass;
		if(verbose && !onePass) printf("Candidate %s failed non targets test\n",primer->seq);
	}
	if (passTests == 8) {
		if(strlen(seqRank)>0){
			rankPrimer(nonDegCandidate , lenSetCand, primer);
		} else {
			primer->tempRanking = 0;
			primer->ranking = 0;
		}
		printPrimer(primer);
		printf("\n");
		fflush(stdout);
		stats->counts [STAT_ALLTESTS] +=1;
		return 1;
	}
	return 0;
}
void checkLoadedPrimers(int forward, Statistics stats) {
	for(int i=0;i< nLoadedPrimers;i++){
		Primer primer = loadedPrimers[i];
		if(primer->forward == forward && testCandidate (primer,stats,0)) {	
			primers[nPrimers] = primer;
			nPrimers++;
		} 
	}
}
void initPrimer(Primer p){
	for(int i=0;i<p->nPositions;i++) {
		p->positions[i]=0;
		p->tempsTargets[i]=-274;
		p->statusTargets[i]=0;
	}
	for(int i=0;i<p->nNonTargets;i++) {
		p->tempsNonTargets[i]=-274;
		p->statusNonTargets[i]=0;
	}
}	
void buildUserPrimer(int forward,Statistics stats) {
	Primer primer = (Primer)malloc (sizeof(TPrimer));
	primer -> id = 1;
	if(forward) {
		strcpy(primer->seq,forwardPrimer);
	} else {
		strcpy(primer->seq,reversePrimer);
	}
	primer -> origin = 'U';
	primer -> forward = forward;
	primer -> nPositions = nTargets; 
	primer -> nNonTargets = nNonTargets; 
	primer -> seqOrig = -1;
	initPrimer(primer);
	if(testCandidate (primer,stats,1)) {	
		primers[nPrimers] = primer;
		nPrimers++;
		if(nPrimers==CANDIDATES_MAX_SIZE) {
			printf("Primers storage capacity reached\n");
		}
	} else {
		free(primer);
	}
	return;

}
void calculatePositionsRange(int nseq, int forward,int * begPos,int * endPos) {
	if(forward) {
		*begPos = beginPosForward;
		*endPos = endPosForward;
	} else {
		*endPos = calculatePosition(nseq,beginPosReverse,0);
		*begPos = calculatePosition(nseq,endPosReverse,0);
	}
	if(*begPos < minLength-1) {
		*begPos = minLength-1;
	}
	if(*endPos > nseq-1) {
		*endPos = nseq-1;
	}	
}
void buildPrimers (int seqCand, int forward, Statistics stats) {
	char * sequence = targets[seqCand];
	int nseq = strlen(sequence);
	int begPos,endPos;
	calculatePositionsRange(nseq,forward,&begPos,&endPos);
	int primerId =1;
	for(int i=begPos;i<=endPos;i++) {
		int ml = maxLength;
		if(ml > i+1 ) {
			ml = i+1;
		}
		for(int length = ml; length>=minLength;length--) {
			Primer primer = (Primer)malloc (sizeof(TPrimer));
			if(buildCandidate(primer->seq,sequence,i,length)) {
				primer -> id = primerId;
				primer -> origin = 'S';
				primer->forward = forward;
				primer -> nPositions = nTargets; 
				primer -> nNonTargets = nNonTargets; 
				primer->seqOrig = seqCand;
				initPrimer(primer);
				if(testCandidate (primer,stats,0)) {	
					primers[nPrimers] = primer;
					primerId++;
					nPrimers++;
					if(nPrimers==CANDIDATES_MAX_SIZE) {
						printf("Primers storage capacity reached\n");
						return;
					}
					float perTC = calculatePercentageTargetsCovered(primer);
					float perNTC = calculatePercentageNonTargetsCovered(primer);
					if(perTC == 100 && perNTC == 0) {
						break;
					}
				} else {
					free(primer);
				}
			}
		}
		int pos = calculatePosition(nseq,i,forward);
		if(pos%200 == 0) printf("Position %d reached\n",pos);	
		fflush(stdout);
	}
}


void runProcess (int forward, Statistics stats) {
	long l2 = time(NULL);
	initialize();
	buildTables ();
	//printSeedsInfo();
	//printNonTargetsSeedsInfo ();
	long l3 = time(NULL);
	printf("Build tables time: %ld\n",(l3-l2));
	fflush(stdout);
	if(nLoadedPrimers>0) {
		checkLoadedPrimers(forward, stats); 
	} else if((forward && strlen(forwardPrimer)>0) || (!forward && strlen(reversePrimer)>0)) {
		buildUserPrimer(forward,stats); 
	} else {
		for(int i=0;i<nTargets && i<nseqForCandidates;i++) {
			buildPrimers(i,forward,stats);
		}
	}
	long l4 = time(NULL);
	printf("Generate primers time: %ld\n",(l4-l3));
	fflush(stdout);
}
int generatePairsForPrimer(Primer p1) {
	char nonDegPrimer [MAX_DEGENERACY][PRIMERS_MAX_SIZE];
	int lenSetPrim = generateNonDegSet (nonDegPrimer,p1->seq);
	for (int i=0;i<nPrimers;i++) {
		Primer p2 = primers[i];
		if (!p2 -> forward ) {
			int j;
			int pass = 1;
			//Same origin test
			if(p1->nPositions != p2->nPositions) {
				pass=0;
			}
			if(p1->nNonTargets != p2->nNonTargets) {
				pass=0;
			}

			//Product length test
			float minPair = 10000;
			float maxPair = -274;
			for(j=0;j<p1->nPositions && pass;j++) {
				if(p1->statusTargets[j] && p2->statusTargets[j]) {
					int productLength = p2->positions[j] + strlen(p2->seq);
					productLength -= (p1->positions[j] - strlen(p1->seq));
					if (productLength < productMinLength || productLength > productMaxLength) {
						pass = 0;
					}
					if(p1->tempsTargets[j] < minPair) {
						minPair = p1->tempsTargets[j];
					}
					if(p2->tempsTargets[j] < minPair) {
						minPair = p2->tempsTargets[j];
					}
					if(p1->tempsTargets[j] > maxPair) {
						maxPair = p1->tempsTargets[j];
					}
					if(p2->tempsTargets[j] > maxPair) {
						maxPair = p2->tempsTargets[j];
					}
				}
			}
			//Same non target coverage test
			for(j=0;j<p1->nNonTargets && pass;j++) {
				if(p1->statusNonTargets[j] && p2->statusNonTargets[j]) {
					pass = 0;
				}
			}	
			//Melting temperature difference test
			if (pass) {
				if(maxPair - minPair > maxPairTempDiff) {
					pass = 0;
				}
			}
			//Melting temperature test to avoid cross hybridization
			if (pass && testComplementarity(nonDegPrimer,lenSetPrim,p2->seq)) {
				pairs [nPairs] = (Pair)malloc(sizeof(TPair));
				pairs [nPairs] -> forwardPrimer = p1;
				pairs [nPairs] -> reversePrimer = p2;
				nPairs++;
				if(nPairs == CANDIDATES_MAX_SIZE) {
					printf("Pairs storage capacity reached\n");
					return 0;
				}
			}
		}
	}
	return 1;
}

void generatePairs () {
	for (int i=0;i<nPrimers;i++) {
		Primer p1 = primers[i];
		if (p1 -> forward) {
			if(!generatePairsForPrimer(p1)) {
				break;
			}
		}		
	}
}

int getPairScore(Pair p) {
	int score = 0;
	for(int k=0;k<p->forwardPrimer->nPositions;k++) {
		if(coveredTargets[k]==0 && p->forwardPrimer->statusTargets[k] && p->reversePrimer->statusTargets[k]) {
			score++;
		}
	}
	return score;
}

int orderPairsSetCover() {
	int i;
	nCoveredTargets = 0;
	for(i=0; i<nPairs && nCoveredTargets<pairs[i]->forwardPrimer->nPositions;i++) {
		int maxScore = getPairScore( pairs[i]);
		int maxPos = i;
		for(int j=i+1; j<nPairs;j++) {
			int nextScore = getPairScore( pairs[j]);
			if(maxScore<nextScore) {
				maxScore = nextScore;
				maxPos = j;
			}
		}
		//Update targets covered
		Pair p = pairs[maxPos];
		for(int k=0;k<pairs[i]->forwardPrimer->nPositions;k++) {
			if(coveredTargets[k]==0 && p->forwardPrimer->statusTargets[k] && p->reversePrimer->statusTargets[k]) {
				coveredTargets[k]=1;
				nCoveredTargets++;
			}
		}
		if(maxPos!=i) {
			Pair p = pairs[i];
			pairs[i] = pairs[maxPos];
			pairs[maxPos] = p;
		}
	}
	return i;
}

void storePrimers(char fileName[]) {
	FILE * file = fopen (fileName, "wt");
	if (file == NULL) {
		printf("File %d could not be created\n",fileName);
		return;
	}
	for(int i=0;i<nPrimers;i++) {
		Primer p = primers[i];
		if(i==0) {
			fprintf(file,"%d,%d,%d\n",p->nPositions,p->nNonTargets,seqPositions);
		}
		fprintf(file,"%d,%s,%d,%d\n",p->id,p->seq,p->forward,p->seqOrig);
		for(int j=0;j<p->nPositions;j++) {
			if(j>0) fprintf(file,",");
			fprintf(file,"%d,%.2f",p->positions[j],p->tempsTargets[j]);
		}
		fprintf(file,"\n");
		for(int j=0;j<p->nNonTargets;j++) {
			if(j>0) fprintf(file,",");
			fprintf(file,"%.2f",p->tempsNonTargets[j]);
		}
		fprintf(file,"\n");
	}
	fclose(file);
}
void loadPrimers(char fileName[]) {
	char line[SEQUENCES_MAX_SIZE];
	FILE * file = fopen (fileName, "rt");
	if (file == NULL) {
		printf("File %d could not be opened\n",fileName);
		return;
	}
	fgets(line,SEQUENCES_MAX_SIZE,file);
	char * tokenptr = strtok(line,",");
	int nTargetsL = atoi(tokenptr);
	tokenptr = strtok(NULL,",");
	int nNonTargetsL = atoi(tokenptr);
	tokenptr = strtok(NULL,",");
	seqPositions  = atoi(tokenptr);
	fgets(line,SEQUENCES_MAX_SIZE,file);
	while(!feof(file)) {
		Primer p = (Primer)malloc (sizeof(TPrimer));
		p->nPositions = nTargetsL;
		p->nNonTargets = nNonTargetsL;
		p->origin = 'L';
		tokenptr = strtok(line,",");
		p->id = atoi(tokenptr);
		tokenptr = strtok(NULL,",");
		strcpy(p->seq,tokenptr);
		tokenptr = strtok(NULL,",");
		p->forward = atoi(tokenptr);
		tokenptr = strtok(NULL,",");
		p->seqOrig = atoi(tokenptr);
		initPrimer(p);
		
		fgets(line,SEQUENCES_MAX_SIZE,file);
		tokenptr = strtok(line,",");
		p -> minTempTargets = 10000;
		p -> maxTempTargets = -274;	
		for(int j=0;j<p->nPositions;j++) {
			p->positions[j] = atoi(tokenptr);
			tokenptr = strtok(NULL,",");
			p->tempsTargets[j] = atof(tokenptr);
			updateTemperatureData(p,p->tempsTargets[j],1);
			tokenptr = strtok(NULL,",");
		}
		fgets(line,SEQUENCES_MAX_SIZE,file);
		tokenptr = strtok(line,",");
		p -> minTempNonTargets = 10000;
		p -> maxTempNonTargets = -274;	
		for(int j=0;j<p->nNonTargets;j++) {
			p->tempsNonTargets[j] = atof(tokenptr);
			updateTemperatureData(p,p->tempsNonTargets[j],0);
			tokenptr = strtok(NULL,",");
		}
		updateAverageTemperature(p);
		loadedPrimers[nLoadedPrimers] = p;
		nLoadedPrimers++; 
		fgets(line,SEQUENCES_MAX_SIZE,file);
	}
	fclose(file);
}	
int validateParameters() {
	//Validate salt method
	if(salMethod != SALT_METHOD_SANTALUCIA && salMethod != SALT_METHOD_OWCZARZY)
	{
		printf("Invalid salt correction method: %d\n", salMethod);
		return 0;
	}
	//Validate primer length
	if (minLength > PRIMERS_MAX_SIZE-1 || maxLength > PRIMERS_MAX_SIZE-1) {
		printf("Max primer length allowed: %d\n", PRIMERS_MAX_SIZE-1);
		return 0;
	}
	//Validate primer positions
	if (beginPosForward < 0 || endPosForward < beginPosForward) {
		printf("Invalid positions range for forward primers ");
	} 
	if (beginPosReverse < 0 || endPosReverse < beginPosReverse) {
		printf("Invalid positions range for reverse primers ");
	} 
	//Validate GCClamp against primer length
	if(gcClamp > minLength) {
		printf("GC Clamp must be less than min primer length ");
	}	
	//Validate tmask
	int count =0;
	for(int i=0; i<strlen(maskTargets);i++) {
		if(maskTargets[i] == '1') {
			count++;
		}
	}
	if(count>SEEDS_MAX_SIZE) {
		printf("Max number of 1 in targets mask allowed: %d\n", SEEDS_MAX_SIZE);
		return 0;
	}
	//Validate nmask
	count =0;
	for(int i=0; i<strlen(maskNonTargets);i++) {
		if(maskNonTargets[i] == '1') {
			count++;
		}
	}
	if(count>SEEDS_MAX_SIZE) {
		printf("Max number of 1 in non targets mask allowed: %d\n", SEEDS_MAX_SIZE);
		return 0;
	}
	//Validate degeneracy
	count =0;
	for(int i=0; i<strlen(degeneracy);i++) {
		if(degeneracy[i] == '4') {
			count++;
		}
	}
	if(pow((float)4,(float)count)>MAX_DEGENERACY) {
		printf("Max degeneracy allowed: %d\n", MAX_DEGENERACY);
		return 0;
	}
	//Validate degeneracy against targets mask
	for(int i=0; i<strlen(degeneracy)&&i<strlen(maskTargets);i++) {
		if(degeneracy[i] == '4' && maskTargets[i] == '1' ) {
			printf("Degeneracy not allowed for required target matching positions\n");
			return 0;
		}
	}

	return 1;
}

int main ( int argc, char ** argv) {
	int sequences=NON_TARGETS_MAX_SIZE;
    	//read command line parameters
	for(int i = 1; i < argc; i++) 
  	{
		if(!strncmp(argv[i], "-tf", strlen("-tf"))) {
			strcpy(targetsFile,argv[++i]);
		} else if(!strncmp(argv[i], "-nf", strlen("-nf"))) {
			strcpy(nonTargetsFile,argv[++i]);
		} else if(!strncmp(argv[i], "-pf", strlen("-pf"))) {
			strcpy(primersFile,argv[++i]);
		} else if(!strncmp(argv[i], "-pof", strlen("-pof"))) {
			strcpy(primersOutputFile,argv[++i]);
		} else if(!strncmp(argv[i], "-minPrimerLength", strlen("-minPrimerLength"))) {
			minLength = atol(argv[++i]);
		} else if(!strncmp(argv[i], "-maxPrimerLength", strlen("-maxPrimerLength"))) {
			maxLength = atol(argv[++i]);
		} else if(!strncmp(argv[i], "-minProdLength", strlen("-minProdLength"))) {
			productMinLength = atol(argv[++i]);
		} else if(!strncmp(argv[i], "-maxProdLength", strlen("-maxProdLength"))) {
			productMaxLength = atol(argv[++i]);
		} else if(!strncmp(argv[i], "-forwardPrimer", strlen("-forwardPrimer"))) {
      			strcpy(forwardPrimer,argv[++i]);
		} else if(!strncmp(argv[i], "-reversePrimer", strlen("-reversePrimer"))) {
      			strcpy(reversePrimer,argv[++i]);
		} else if(!strncmp(argv[i], "-beginPosForward", strlen("-beginPosForward"))) {
      			beginPosForward = atol(argv[++i]);
    		} else if(!strncmp(argv[i], "-endPosForward", strlen("-endPosForward"))) {
      			endPosForward = atol(argv[++i]);
		} else if(!strncmp(argv[i], "-beginPosReverse", strlen("-beginPosReverse"))) {
      			beginPosReverse = atol(argv[++i]);
    		} else if(!strncmp(argv[i], "-endPosReverse", strlen("-endPosReverse"))) {
      			endPosReverse = atol(argv[++i]);
    		} else if(!strncmp(argv[i], "-tmask", strlen("-tmask"))) {
      			strcpy(maskTargets,argv[++i]);
    		} else if(!strncmp(argv[i], "-nmask", strlen("-nmask"))) {
      			strcpy(maskNonTargets,argv[++i]);
    		} else if(!strncmp(argv[i], "-dmask", strlen("-dmask"))) {
      			strcpy(degeneracy,argv[++i]);
    		} else if(!strncmp(argv[i], "-sequences", strlen("-sequences"))) {
      			sequences = atol(argv[++i]);
    		} else if(!strncmp(argv[i], "-numSourceSeq", strlen("-numSourceSeq"))) {
      			nseqForCandidates = atol(argv[++i]);
    		} else if(!strncmp(argv[i], "-minCoverageTargets", strlen("-minCoverageTargets"))) {
      			minCoverageTargets = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-maxCoverageNonTargets", strlen("-maxCoverageNonTargets"))) {
      			maxCoverageNonTargets = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-maxSelfScore", strlen("-maxSelfScore"))) {
      			maxSelfScore = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-maxEndScore", strlen("-maxEndScore"))) {
      			maxEndScore = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-minGCContent", strlen("-minGCContent"))) {
      			minGCContent = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-maxGCContent", strlen("-maxGCContent"))) {
      			maxGCContent = atof(argv[++i]);
		} else if(!strncmp(argv[i], "-gcClamp", strlen("-gcClamp"))) {
      			gcClamp = atoi(argv[++i]);
		} else if(!strncmp(argv[i], "-maxPolyX", strlen("-maxPolyX"))) {
      			maxPolyX = atoi(argv[++i]);
		} else if(!strncmp(argv[i], "-primersConc", strlen("-primersConc"))) {
      			concPrimers = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-templateConc", strlen("-templateConc"))) {
      			concSequences = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-saltConc", strlen("-saltConc"))) {
      			salt = atof(argv[++i]);
		} else if(!strncmp(argv[i], "-saltCorrMethod", strlen("-saltCorrMethod"))) {
      			salMethod = atoi(argv[++i]);
		} else if(!strncmp(argv[i], "-minTempTargets", strlen("-minTempTargets"))) {
      			minTempTargets = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-maxTempTargets", strlen("-maxTempTargets"))) {
      			maxTempTargets = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-maxTempNonTargets", strlen("-maxTempNonTargets"))) {
      			maxTempNonTargets = atof(argv[++i]);
		} else if(!strncmp(argv[i], "-deltaTempNonTargets", strlen("-deltaTempNonTargets"))) {
      			deltaTempNonTargets = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-maxPairTempDiff", strlen("-maxPairTempDiff"))) {
      			maxPairTempDiff = atof(argv[++i]);
    		} else if(!strncmp(argv[i], "-primersLabel", strlen("-primersLabel"))) {
      			strcpy(primersLabel,argv[++i]);
    		} else if(!strncmp(argv[i], "-full_stats", strlen("-full_stats"))) {
      			fullStats = 1;
    		} else if(!strncmp(argv[i], "-seqRank", strlen("-seqRank"))) {
      			strcpy(seqRank,argv[++i]);
    		} else /* unrecognized parameter */ {
      			print_usage( argv[0] );
      			return  -1;
    		}
	}
	printParameters(  );
	//Preprocess parameters
	if(!strcmp(maskNonTargets,"NONE")) {
		maskNonTargets[0]=0;
	}
	if(!strcmp(forwardPrimer,"NONE")) {
		forwardPrimer[0]=0;
	}
	if(!strcmp(reversePrimer,"NONE")) {
		reversePrimer[0]=0;
	}
	if(!strcmp(seqRank,"NONE")) {
		seqRank[0]=0;
	}

	depureSequence(forwardPrimer);
	depureSequence(reversePrimer);
	depureSequence(seqRank);
 
	if (!validateParameters()) {
		return -1;
	}

	long l1 = time(NULL);
	paramsPrimerSeq = new CNNParams();
	paramsPrimerSeq->InitParams(concPrimers,concSequences,salt,salMethod);

	set_dpal_args(&local_args);
	local_args.flag = DPAL_LOCAL;               
	dpal_set_ambiguity_code_matrix(&local_args);

	set_dpal_args(&end_args);
	end_args.flag = DPAL_GLOBAL_END;
	dpal_set_ambiguity_code_matrix(&end_args);

	loadFiles(sequences);
	if(strlen(primersFile)>0) {
		loadPrimers(primersFile);
	} else {
		seqPositions = getSequenceMaxLength(targets,nTargets); 
	}
	//printTargets();
	printf("Targets: %d\n",nTargets);
	printf("Non Targets: %d\n",nNonTargets);
	long l2 = time(NULL);
	printf("Load time: %ld\n",(l2-l1));
	Statistics stats = (Statistics) malloc(sizeof (TStatistics));
	for(int i = 0; i < STAT_COUNTS; i++) {
		stats -> counts [i] = 0;
	}
	for(int i = 0; i < TARGETS_MAX_SIZE; i++) {
		stats -> failTargetsCounts[i] = 0;
		
	}

	printf("Running process for forward primers\n");
	fflush(stdout);
	runProcess(1,stats);
	reverseComplementSequences(targets,nTargets);
	reverseComplementSequences(nontargets,nNonTargets);
	char tmp [SEQUENCES_MAX_SIZE];
	calculateReverseComplement(tmp,seqRank);
	strcpy(seqRank,tmp);
	printf("Running process for reverse primers\n");
	runProcess(0,stats);
	if(strlen(primersOutputFile)>0) {
		storePrimers(primersOutputFile);
	}
	//printPrimers();
	generatePairs();
	if (nPairs > 0 && minCoverageTargets<100) {
		for(int i = 0; i < TARGETS_MAX_SIZE; i++) {
			coveredTargets [i] = 0;
		}
		nCoveredTargets = 0;
		int pairsFullCoverage = orderPairsSetCover();
		if(nCoveredTargets == pairs[0]->forwardPrimer->nPositions) {
			printf("First %d pairs cover all targets\n",pairsFullCoverage);
		} else {
			printf("%d targets non covered by any pair\n",pairs[0]->forwardPrimer->nPositions-nCoveredTargets);
		}
	}
	printPairs();
	printStatistics(stats);	
}
