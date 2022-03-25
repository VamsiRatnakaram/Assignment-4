#ifndef _WIRE_ROUTE_H
#define _WIRE_ROUTE_H

typedef struct bend{
    int x;
    int y;
} bend_t;

typedef struct totalCost {
    int maxValue;
    int cost;
} total_cost_t;

typedef struct { /* Define the data structure for wire here */
    int start_x;
    int start_y;
    int end_x;
    int end_y;
    int numBends;
    int bend0x;
    int bend0y;
    int bend1x;
    int bend1y;
} wire_t;

typedef int cost_t;

// Perform computation, including reading/writing output files
int compute(int procID, int nproc, char *inputFilename, double prob, int numIterations);

// Read input file
void readInput(char *inputFilename);

// Write cost array file based on input filename
void writeCost(char *inputFilename, int nproc);

// Write wire output file based on input filename
void writeOutput(char *inputFilename, int nproc);

#endif // _WIRE_ROUTE_H
