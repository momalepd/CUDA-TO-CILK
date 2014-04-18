/** Single source shortest paths -*- C++ -*-
	- Modified Bellman Ford algorithm
	- Implementation of CUDA version of GPU to Cilk Plus on Multicore Processor
	- Lonestargpu benchmark is used to migrate from CUDA to Cilk Plus
	- This is a plain transformation of CUDA BFS to Cilk Plus.
*/

#include "_lonestargpu.h"
#include "_sssp_worklist.h"

//-------------------------------------------------------------------------------------------------
// - Veryfing the final distances

void dverifysolution(foru *dist, Graph graph, cilk::reducer_opadd<int> &nerr) {
	cilk_for(int nn = 0 ;  nn < graph.nnodes ; nn++) {
		unsigned int nsrcedges = graph.getOutDegree(nn);
		for (unsigned ii = 0; ii < nsrcedges; ++ii) {
			unsigned int u = nn;
			unsigned int v = graph.getDestination(u, ii);
			foru wt = graph.getWeight(u, ii);
			if (wt > 0 && dist[u] + wt < dist[v]) {
				++*nerr;
			}
		}
	}	
}

//-------------------------------------------------------------------------------------------------
// - Writing final output to the file

void print_output(const char *filename, foru *hdist, foru *dist, Graph graph) {
	CUDA_SAFE_CALL(cudaMemcpy(hdist, dist, graph.nnodes * sizeof(foru), cudaMemcpyDeviceToHost));

  	printf("Writing output to %s\n", filename);
  	FILE *o = fopen(filename, "w");

  	for(int i = 0; i < graph.nnodes; i++) {
    	fprintf(o, "%d: %d\n", i, hdist[i]);
  	}

  	fclose(o);
}

//-------------------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
	foru *dist;						// Contains final shortest distance from src node
	Graph hgraph;					// Graph data structure to parse graph file
	cilk::reducer_opadd<int> nerr;  // to store total err after the BFS for veryfication

	if (argc != 2) {
		printf("Usage: %s <graph>\n", argv[0]);
		exit(1);
	}

	hgraph.read(argv[1]);  // Read give file into graph Data Structure
	//hgraph.optimize();   // Optimizing graph

	dist = (foru *)malloc(graph.nnodes * sizeof(foru));  

	// Main sssp function
	sssp(dist, hgraph);

	// Veryfing the final distances
	printf("verifying.\n");
	verifysolution (dist, graph, nerr);
	printf("\tno of errors = %d.\n", nerr.get_value());
	
	print_output("sssp-output.txt", hdist, dist, graph);
	// cleanup left to the OS.

	return 0;
}
