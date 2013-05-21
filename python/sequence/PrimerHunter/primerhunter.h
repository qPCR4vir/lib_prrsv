#ifndef _PRIMERHUNTER_H
#define _PRIMERHUNTER_H

#define TARGETS_MAX_SIZE 2000
#define NON_TARGETS_MAX_SIZE 10000
#define SEQUENCES_MAX_SIZE 100000
#define PRIMERS_MAX_SIZE 100
#define CANDIDATES_MAX_SIZE 100000
#define SEEDS_MAX_SIZE 9
#define SEEDSTABLE_TAM 1000000
#define MAX_DEGENERACY 256

#define STAT_COUNTS 12
#define STAT_TOTAL_CANDIDATES 0
#define STAT_TEST_DEGENERACY 1
#define STAT_TEST_GCCONTENT 2
#define STAT_TEST_GCCLAMP 3
#define STAT_TEST_MAX_POLY_X 4
#define STAT_TEST_SELFCOMP 5
#define STAT_TEST_TMSELF 6
#define STAT_TEST_TARGETS 7
#define STAT_TEST_NONTARGETS 8
#define STAT_ALLTESTS 9
#define STAT_FAILED_SELFCOMP_LOW 10
#define STAT_FAILED_SELFCOMP_HIGH 11

#define DEF_MIN_PRIMER_LENGTH 20
#define DEF_MAX_PRIMER_LENGTH 25
#define DEF_MIN_PROD_LENGTH 75
#define DEF_MAX_PROD_LENGTH 200
#define DEF_FORWARD_PRIMER "NONE"
#define DEF_REVERSE_PRIMER "NONE"
#define DEF_BEGINPOS_FORWARD 0
#define DEF_ENDPOS_FORWARD SEQUENCES_MAX_SIZE
#define DEF_BEGINPOS_REVERSE 0
#define DEF_ENDPOS_REVERSE SEQUENCES_MAX_SIZE
#define DEF_MASK_TARGETS "11"
#define DEF_MASK_NONTARGETS "NONE"
#define DEF_DEGENERACY "1"
#define DEF_SEQ_FOR_CANDS 1
#define DEF_MIN_COVERAGE_TARGETS 100
#define DEF_MAX_COVERAGE_NONTARGETS 0
#define DEF_MAX_SELF_SCORE 800
#define DEF_MAX_END_SCORE 300
#define DEF_MIN_GCCONTENT 25
#define DEF_MAX_GCCONTENT 75
#define DEF_GCCLAMP 0
#define DEF_MAX_POLY_X 5
#define DEF_CONC_PRIMERS 0.0000008
#define DEF_CONC_SEQUENCES 0
#define DEF_SALT 0.05
#define DEF_MIN_TEMP_TARGETS 40
#define DEF_MAX_TEMP_TARGETS 70
#define DEF_MAX_TEMP_NONTARGETS 40
#define DEF_DELTA_TEMP_NONTARGETS 0
#define DEF_MAX_PAIR_TEMP_DIFF 40
#define DEF_PRIMERS_LABEL "P"
#define DEF_FULL_STATS 0
#define DEF_SEQ_RANK "NONE"

typedef struct occurance {
	int sequenceId; 
	int sequencePosition; //Position of the last character
	occurance * next;
} TOccurance, *Occurance;

typedef struct {
	Occurance head;
	Occurance tail;
	int size;
} TOccuranceList, *OccuranceList;

typedef struct {
	char seed [SEEDS_MAX_SIZE];
        OccuranceList occurancesList;
} TSeed, * Seed;

typedef struct {
	int id;
	char origin; //Tells where the primer comes from. S for sequence, U for user, L for loaded  
	char seq [PRIMERS_MAX_SIZE]; //Primer sequence in 5' - 3' orientation
	int forward; //1 for forward. 0 for reverse
	int positions [TARGETS_MAX_SIZE];
	float tempsTargets [TARGETS_MAX_SIZE];
	int statusTargets[TARGETS_MAX_SIZE]; //1 if the primer hybridizes the target, 0 otherwise 
	int nPositions;
	float tempsNonTargets [NON_TARGETS_MAX_SIZE];
	int statusNonTargets[NON_TARGETS_MAX_SIZE];//1 if the primer hybridizes the non target, 0 otherwise
	int nNonTargets;
	int seqOrig;
	float minTempTargets;
	float avgTempTargets;
	float maxTempTargets;
	float minTempNonTargets;
	float maxTempNonTargets;
	float tempRanking;
	float ranking; //Proportion of targets with lower melting temp than the temperature between the primer and the ranking sequence 
} TPrimer, * Primer;

typedef struct {
	Primer forwardPrimer;
	Primer reversePrimer;
} TPair, * Pair;

typedef struct {
	int counts [STAT_COUNTS];
	int failTargetsCounts [TARGETS_MAX_SIZE];
} TStatistics, *Statistics;

#endif
