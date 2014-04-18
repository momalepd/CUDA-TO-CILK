/* Breadth First Search
	- Implementation of CUDA version of GPU to Cilk Plus on Multicore Processor
	- Lonestargpu benchmark is used to migrate from CUDA to Cilk Plus
	- This is a plain transformation of CUDA BFS to Cilk Plus.
*/


#include "lonestargpu.h"
#include "bfs_worklistc.h"

//--------------------------------------------------------------------------------------------
// -Verifying the solution

void verifysolution(foru *dist, Graph graph, cilk::reducer_opadd<int> &nerr) {
	
	cilk_for(int nn = 0 ;  nn < graph.nnodes ; nn++) {
		unsigned int nsrcedges = graph.getOutDegree(nn);
		for (unsigned ii = 0; ii < nsrcedges; ++ii) {
			unsigned int u = nn;
			unsigned int v = graph.getDestination(u, ii);
			foru wt = 1;

			if (wt > 0 && dist[u] + wt < dist[v]) {
			  //printf("%d %d %d %d\n", u, v, dist[u], dist[v]);
			  nerr += 1;
			}
		}
	  }	
}




//--------------------------------------------------------------------------------------------
// - Writing the solution to the output file "bfs-output.txt"


void write_solution(const char *fname, Graph &graph, foru *dist)
{
  
  printf("Writing solution to %s\n", fname);
  FILE *f = fopen(fname, "w");
  // formatted like Merrill's code for comparison
  fprintf(f, "Computed solution (source dist): [");

  for(int node = 0; node < graph.nnodes; node++)
    {
      fprintf(f, "%d:%d\n ", node, dist[node]);
    }

  fprintf(f, "]");

  free(dist);
}

//--------------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
	
	unsigned intzero = 0;
	Graph hgraph;    // Graph data structure to parse graph file
	foru *dist;		 // contains result (level no. of each vertex after the BFS)
	cilk::reducer_opadd<int> nerr;  // to store total err after the BFS for veryfication


	if (argc != 2) {
		printf("Usage: %s <graph>\n", argv[0]);
		exit(1);
	}

	hgraph.read(argv[1]);  // Read give file into graph Data Structure

	dist = (foru *)malloc(graph.nnodes * sizeof(foru));

	bfs(hgraph, dist);   // Main BFS execution

	printf("verifying.\n");
	verifysolution (dist, graph, nerr);
	printf("\tno of errors = %d.\n", nerr.get_value());

	write_solution("bfs-output.txt", hgraph, dist);

	// cleanup left to the OS.

	return 0;
}


