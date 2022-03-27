#include "wireroute.h"
#include "mpi.h"
#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libgen.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <linux/limits.h>
#include <string>

using namespace std;

static void update_route(wire_t wire,int *costs,int dim_x,int dim_y,int updateVal){
        int start_x = wire.start_x;
        int start_y = wire.start_y;
        int end_x,end_y;
        for(int i=0; i < wire.numBends+1; i++){
            if(i==wire.numBends){
                end_x = wire.end_x;
                end_y = wire.end_y;
            }else{
                end_x = (i != 1) ? wire.bend0x : wire.bend1x;
                end_y = (i != 1) ? wire.bend0y : wire.bend1y;
            }
            if(start_x==end_x){
                for(int j = min(start_y,end_y); j < max(end_y,start_y)+1; j++){
                    int temp_bendy = (i == 1) ? wire.bend0y : wire.bend1y; 
                    if (i>0 && temp_bendy==j) {
                        continue;
                    }
                    costs[j*dim_y+start_x]+=updateVal;
                } 
                start_y=end_y;
            }else{
                for(int k = min(start_x,end_x); k < max(end_x,start_x)+1; k++){
                    int temp_bendx = (i == 1) ? wire.bend0x : wire.bend1x; 
                    if (i>0 && temp_bendx==k) {
                        continue;
                    }
                    costs[start_y*dim_y+k]+=updateVal;
                }
                start_x=end_x;
            }
        }
}

static void initialize(wire_t *wires, int *costs, int dim_x,int dim_y,int num_wires) {
    //looks messy maybe optimize
    for (int i = 0; i < num_wires; i++){
        wire_t *init_wire = &wires[i];
        if(init_wire->start_x==init_wire->end_x || init_wire->start_y==init_wire->end_y){
            init_wire->numBends=0;
        }else{
            init_wire->numBends=1;
            init_wire->bend0x = init_wire->end_x;
            init_wire->bend0y = init_wire->start_y;
        }
        update_route(*init_wire,costs,dim_x,dim_y,1);
    }
}

static void createCostMap(wire_t *old_wires,wire_t *wires, int *costs, int dim_x,int dim_y,int num_wires) {
    for (int i = 0; i < num_wires; i++){
        update_route(old_wires[i],costs,dim_x,dim_y,-1);
        update_route(wires[i],costs,dim_x,dim_y,1);
    }
}

static total_cost_t calculateCost(wire_t curr, int *costs, int dim_x, int dim_y) {
    int maxVal=0;
    int start_x=curr.start_x;
    int start_y=curr.start_y;
    int end_x,end_y;
    int currCost = 0;
    int numBends = curr.numBends;
    int x=dim_x;
    int y=dim_y;
    for(int i=0; i < numBends+1; i++){
        if(i==numBends){
            end_x=curr.end_x;
            end_y=curr.end_y;
        }else{
            end_x = (i == 0) ? curr.bend0x : curr.bend1x;
            end_y = (i == 0) ? curr.bend0y : curr.bend1y;
        }
        if(start_x==end_x){
            int start = min(start_y,end_y);
            int end=max(end_y,start_y);
            int index=start*y+start_x;
            int end_index=end*y+start_x;
            int j=index;
            for(; j < end_index+1; j+=2*y){
                int temp_bendy = (i == 1) ? curr.bend0y : curr.bend1y; // TO-DO check correctness of this line
                if (!(i && temp_bendy==j)) {
                    int c1=costs[j];
                    int c2=costs[j+y];
                    currCost+=c1+c2+1;
                    maxVal=max(maxVal,max(c1,c2));
                }
            } if(j!=end_index){
                int c=costs[end_index];
                currCost+=c;
                maxVal=max(maxVal,c);
            }
            start_y=end_y;
        }else{
            int start=min(start_x,end_x);
            int end=max(end_x,start_x);
            int index= start_y*y+start;
            int end_index=start_y*y+end;
            int k=index;
            for(; k < end_index+1; k+=2){
                int temp_bendx = (i == 1) ? curr.bend0x : curr.bend1x; // TO-DO check correctness of this line
                if (!(i && temp_bendx==k)) {
                    int c1=costs[k];
                    int c2=costs[k+1];
                    currCost+=c1+c2+1;
                    maxVal=max(maxVal,max(c1,c2));
                }
            }if(k!=end_index){
                int c=costs[end_index];
                currCost+=c;
                maxVal=max(maxVal,c);
            }
            start_x=end_x;
        }
    }
    total_cost_t total_cost;
    total_cost.cost=currCost;
    total_cost.maxValue=maxVal;
    return total_cost;
}

void defineWireStruct(MPI_Datatype *tstype) {
    const int count = 9;
    int          blocklens[count] = {1,1,1,1,1,1,1,1,1};
    MPI_Datatype types[9] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Aint     disps[count] = {offsetof(wire_t,start_x), offsetof(wire_t,start_y), offsetof(wire_t,end_x), offsetof(wire_t,end_y), offsetof(wire_t,numBends), offsetof(wire_t,bend0x), offsetof(wire_t,bend0y), offsetof(wire_t,bend1x), offsetof(wire_t,bend1y)};

    MPI_Type_create_struct(count, blocklens, disps, types, tstype);
    MPI_Type_commit(tstype);
}

void defineWireValidStruct(MPI_Datatype *tstype) {
    const int count = 11;
    int          blocklens[count] = {1,1,1,1,1,1,1,1,1,1,1};
    MPI_Datatype types[11] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Aint     disps[count] = {offsetof(wireValid_t,start_x), offsetof(wireValid_t,start_y), offsetof(wireValid_t,end_x), offsetof(wireValid_t,end_y), offsetof(wireValid_t,numBends), offsetof(wireValid_t,bend0x), offsetof(wireValid_t,bend0y), offsetof(wireValid_t,bend1x), offsetof(wireValid_t,bend1y), offsetof(wireValid_t,valid), offsetof(wireValid_t,wireIndex)};

    MPI_Type_create_struct(count, blocklens, disps, types, tstype);
    MPI_Type_commit(tstype);
}

static wire_t updateHelper(wire_t InWire, int *costs, int dim_x, int dim_y, int random_prob) {
    int numBends = InWire.numBends;
    int start_x = InWire.start_x;
    int start_y = InWire.start_y;
    int end_x = InWire.end_x;
    int end_y = InWire.end_y;
    if (numBends == 0) {
        return InWire;
    }

    int sign_x=1,sign_y=1;
    if(start_x > end_x){
        sign_x=-1;
    }
    if(start_y > end_y){
        sign_y=-1;
    }

    // Create Random Path
    // Horizontal first
    wire_t hWire = InWire;
    hWire.bend0x = (rand() % (abs(end_x-start_x)))*sign_x + start_x;
    hWire.bend0y = start_y;
    if (hWire.bend0x == end_x) {
        hWire.numBends = 1;
    }
    else {
        hWire.numBends = 2;
        hWire.bend1x = hWire.bend0x;
        hWire.bend1y = end_y;
    }
    // Vertical First
    wire_t vWire = InWire;
    vWire.bend0y = (rand() % (abs(end_y-start_y)))*sign_y + start_y;
    vWire.bend0x = start_x;
    if (vWire.bend0y == end_y) {
        vWire.numBends = 1;
    }
    else {
        vWire.numBends = 2;
        vWire.bend1y = vWire.bend0y;
        vWire.bend1x = end_x;
    }
    int h_or_v = rand() % 2;
    int randomProb = rand() % 100;
    // Replace best wire with random wire
    if (randomProb < random_prob) {
        wire_t randomWire = (h_or_v) ? hWire : vWire;
        return randomWire;
    }

    total_cost_t currCost = calculateCost(InWire, costs, dim_x, dim_y);
    wire_t newWire = InWire;
    wire_t bestWire = InWire;
    newWire.numBends = 0;

    // Check All Possible Paths
    for (int j = 0; j < abs(end_x - start_x) + abs(end_y - start_y); j+=1) {
        if (j < abs(end_x - start_x)) { // Horizontal Path
            newWire.bend0x = start_x + sign_x*(j + 1);
            newWire.bend0y = start_y;
            if (newWire.bend0x == end_x) {
                // One Bend Case
                newWire.numBends = 1;
            }
            else {
                // Two Bend Case
                newWire.bend1x = newWire.bend0x;
                newWire.bend1y = end_y;
                newWire.numBends = 2;
            }
        }
        else { // Once done with all horizontal paths, we can compute vertical paths
            int modified_j = j - abs(end_x - start_x);
            newWire.bend0y = start_y + sign_y*(modified_j + 1);
            newWire.bend0x = start_x;
            if ( newWire.bend0y == end_y) {
                // One Bend Case
                newWire.numBends = 1;
            }
            else {
                // Two Bend Case
                newWire.bend1y = newWire.bend0y;
                newWire.bend1x = end_x;
                newWire.numBends = 2;
            }
            // Check if newWire is better than oldWire and replace if so
        }
        total_cost_t newCost = calculateCost(newWire, costs, dim_x, dim_y);
        if(newCost.maxValue < currCost.maxValue){
            currCost = newCost;
            bestWire = newWire; 
        }else if (newCost.maxValue == currCost.maxValue && newCost.cost < currCost.cost) {
            currCost = newCost;
            bestWire = newWire;
        }
    }
    return bestWire;
}

static wire_t wireValidToWire(wireValid_t wireValid) {
    wire_t temp;
    temp.start_x = wireValid.start_x;
    temp.start_y = wireValid.start_y;
    temp.end_x = wireValid.end_x;
    temp.end_y = wireValid.end_y;
    temp.numBends = wireValid.numBends;
    temp.bend0x = wireValid.bend0x;
    temp.bend0y = wireValid.bend0y;
    temp.bend1x = wireValid.bend1x;
    temp.bend1y = wireValid.bend1y;
    return temp;
}

static wireValid_t wireToWireValid(wire_t wire, int valid, int index) {
    wireValid_t temp;
    temp.valid = valid;
    temp.wireIndex = index;
    temp.start_x = wire.start_x;
    temp.start_y = wire.start_y;
    temp.end_x = wire.end_x;
    temp.end_y = wire.end_y;
    temp.numBends = wire.numBends;
    temp.bend0x = wire.bend0x;
    temp.bend0y = wire.bend0y;
    temp.bend1x = wire.bend1x;
    temp.bend1y = wire.bend1y;
    return temp;
}

static void update(wire_t *wires, int *costs, int dim_x, int dim_y, int num_wires, int random_prob, int procID, int nproc, int *displacements, int *counts, wireValid_t *sendBuf, wireValid_t *recvBuf, int batch_size) {
    // Find better cost for each wire
    (void)procID;
    int rem = num_wires%nproc;
    MPI_Datatype wireStruct;
    defineWireStruct(&wireStruct);
    MPI_Datatype wireValidStruct;
    defineWireValidStruct(&wireValidStruct);
    wire_t newWire;
    wireValid_t *newValidWire = (wireValid_t *)calloc(batch_size, sizeof(wireValid_t));
    int f = 0;
    for (int i = displacements[procID]; i < displacements[procID] + counts[procID]; i+=1) {
        // Wires we are modifying
        wire_t oldWire = wires[i];

        update_route(oldWire,costs,dim_x,dim_y,-1);

        newWire = updateHelper(oldWire, costs, dim_x, dim_y, random_prob);

        update_route(newWire,costs,dim_x,dim_y,1);
        wires[i] = newWire;
        newValidWire[f%batch_size] = wireToWireValid(newWire, 1, i);

        // If we are doing remainder wires
        if (procID < rem && i == displacements[procID] + counts[procID] - 1) {
            break;
        }

        if (f%batch_size == 0) {
            // Communication Here
            memcpy(sendBuf, recvBuf, nproc * batch_size * sizeof(wireValid_t));
            memcpy((wireValid_t *)(&sendBuf[procID*batch_size]), newValidWire, batch_size * sizeof(wireValid_t));
            
            // Asynch Send
            MPI_Request request;
            int sendProc = (procID == 0) ? nproc - 1 : procID - 1;
            MPI_Isend(sendBuf, batch_size*nproc, wireValidStruct, sendProc, procID, MPI_COMM_WORLD, &request);
            // Synch Recieve
            int recvProc = (procID == nproc - 1) ? 0 : procID + 1;
            MPI_Recv(recvBuf, batch_size*nproc, wireValidStruct, recvProc, recvProc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Wait(&request,MPI_STATUS_IGNORE);
            // Update Wires that you recieved
            for (int k = 0; k < nproc*batch_size; k++) {
                if (k/batch_size == procID) {
                    continue;
                }
                if (recvBuf[k].valid) {
                    wire_t temp = wireValidToWire(recvBuf[k]);
                    int index = recvBuf[k].wireIndex;
                    update_route(wires[index], costs, dim_x, dim_y,-1);
                    wires[index] = temp;
                    update_route(wires[index], costs, dim_x, dim_y,1);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        f++;
    }
    free(newValidWire);
}

// Perform computation, including reading/writing output files
double compute(int procID, int nproc, char *inputFilename, double prob, int numIterations) {
    // TODO Implement code here
    // TODO Decide which processors should be reading/writing files
    const int root = 0; // Set the rank 0 process as the root process
    int tag = 0;        // Use the same tag
    double startTime;
    double endTime;
    MPI_Status status;
    MPI_Datatype wireStruct;

    FILE *input;

    input = fopen(inputFilename, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", inputFilename);
        return -1;
    }

    int dim_x, dim_y;
    int num_of_wires;

    int i = fscanf(input, "%d %d\n", &dim_x, &dim_y);
    int j = fscanf(input, "%d\n", &num_of_wires);
    (void)i;
    (void)j;
    wire_t *wires = (wire_t *)calloc(num_of_wires, sizeof(wire_t));
    /* Read the grid dimension and wire information from file */

    cost_t *costs = (cost_t *)calloc(dim_x * dim_y, sizeof(cost_t));
    /* Initialize cost matrix */

    // For Asynch Send and Recieve

    // Read Input File and Initialize Arrays
    if (procID == root) {
        /* Initailize additional data structures needed in the algorithm */
        for (int i = 0; i < num_of_wires; i++) {
            int k = fscanf(input, "%d %d %d %d\n", &(wires[i].start_x), &(wires[i].start_y), &(wires[i].end_x), &(wires[i].end_y));
            (void)k;
        }
        /* Conduct initial wire placement */
        initialize(wires,costs,dim_x,dim_y,num_of_wires);
    }
    int batch_size = num_of_wires/nproc;

    wireValid_t *sendBuf = (wireValid_t *)calloc(nproc*batch_size, sizeof(wireValid_t));
    wireValid_t *recvBuf = (wireValid_t *)calloc(nproc*batch_size, sizeof(wireValid_t));


    // Create MPI Structs to send using MPI
    defineWireStruct(&wireStruct);
    // Figuring out scatterv and gatherv distribution
    int *counts = (int*)calloc(nproc, sizeof(int)); // array describing how many elements to send to each process
    int rem = (num_of_wires)%nproc; // elements remaining after division among processes
    int *displacements = (int*)calloc(nproc, sizeof(int)); // array describing the displacements where each segment begins
    int sum = 0;                // Sum of counts. Used to calculate displacements

    // calculate send counts and displacements
    for (int k = 0; k < nproc; k++) {
        counts[k] = num_of_wires/nproc;
        if (rem > 0) {
            counts[k]++;
            rem--;
        }

        displacements[k] = sum;
        sum += counts[k];
    }

    // StartTime after intialization
    startTime = MPI_Wtime();

    // Bcast wires and costs
    MPI_Bcast(wires, num_of_wires, wireStruct, root, MPI_COMM_WORLD);
    MPI_Bcast(costs, dim_x * dim_y, MPI_INT, root, MPI_COMM_WORLD);

    for (int i = 0; i < numIterations; i++) {
        // update function HERE
        update(wires, costs, dim_x, dim_y, num_of_wires, (int)(100*prob), procID, nproc, displacements, counts, sendBuf, recvBuf, batch_size);
    }

    // EndTime before I/O
    endTime = MPI_Wtime();

    // Output Files I/O
    if (procID == root) {
        /* Write wires and costs to files */
        char resolved_path[PATH_MAX];
        char *a;
        a = realpath(inputFilename, resolved_path);
        (void)a;
        char *base = basename(resolved_path);
        std::string baseS = std::string(base);
        size_t lastindex = baseS.find_last_of("."); 
        std::string rawname = baseS.substr(0, lastindex); 

        std::stringstream OutputCosts;
        OutputCosts << "outputs//costs_" << rawname.c_str() << "_" << nproc << ".txt";
        std::string OutputCostsFile = OutputCosts.str();

        std::stringstream OutputWires;
        OutputWires << "outputs//output_" << rawname.c_str() << "_" << nproc << ".txt";
        std::string OutputWiresFile = OutputWires.str();

        const char *ocf = OutputCostsFile.c_str();
        FILE *costFile = fopen(ocf, "w");
        if (!costFile) {
            printf("Unable to open file: %s.\n", ocf);
            return -1;
        }
        const char *owf = OutputWiresFile.c_str();
        FILE *outFile = fopen(owf, "w");
        if (!outFile) {
            printf("sad\n");
            printf("Unable to open file: %s.\n", owf);
            return -1;
        }

        // Write to cost file
        int maxCost = 0;
        int sumOfSquares = 0;
        fprintf(costFile, "%d %d\n", dim_x, dim_y);
        for(int i = 0; i < dim_y; i++){
            for(int j = 0; j < dim_x; j++){
                int temp = costs[i*dim_y+j];
                fprintf(costFile, "%d ", temp);
                if(maxCost < temp) maxCost = temp;
                sumOfSquares += temp*temp;
            }
            fprintf(costFile, "\n");
        }

        // Write to output wire file
        fprintf(outFile, "%d %d\n", dim_x, dim_y);
        fprintf(outFile, "%d\n", num_of_wires);
        for (int i = 0; i < num_of_wires; i++) {
            wire_t curr = wires[i];
            int start_x = curr.start_x;
            int start_y = curr.start_y;
            int end_x,end_y;
            for(int i = 0; i < curr.numBends+1; i++){
                if(i==curr.numBends){
                    end_x=curr.end_x;
                    end_y=curr.end_y;
                }else{
                    end_x = (i != 1) ? curr.bend0x : curr.bend1x;
                    end_y = (i != 1) ? curr.bend0y : curr.bend1y;
                }
                if(start_x==end_x){
                    int sign = start_y < end_y ? 1 : -1;
                    for(int j = 0; j < abs(end_y-start_y)+1; j++){
                        if(i>0 && j==0) continue;
                        fprintf(outFile, "%d %d ", start_x, start_y+j*(sign));
                    }
                    start_y=end_y;
                }else{
                    int sign = start_x < end_x ? 1 : -1;
                    for(int j = 0; j < abs(end_x-start_x)+1; j++){
                        if(i>0 && j==0) continue;
                        fprintf(outFile, "%d %d ", start_x+j*(sign),start_y);
                    }
                    start_x=end_x;

                }
            }
            fprintf(outFile, "\n");
        }

        printf("MaxCost: %d\n", maxCost);
        printf("Sum Of Squares: %d\n", sumOfSquares);

        // Close all files
        free(counts);
        free(displacements);
        free(sendBuf);
        free(recvBuf);
        free(wires);
        free(costs);
        fclose(input);
        fclose(costFile);
        fclose(outFile);
    }
    return endTime-startTime;
}
