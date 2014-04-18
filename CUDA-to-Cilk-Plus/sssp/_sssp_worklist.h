#define MAXDIST		100

#define AVGDEGREE	2.5
#define WORKPERTHREAD	1

unsigned int NVERTICES;

#include "_worklistc.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <list>


#define cilk::reducer_list_append<int> Cilk_Reduce 
#define std::list<int> List



//--------------------------------------------------------------------------------------------------
// - Initializing all the distances to infinity

void initialize(foru *dist, unsigned int nv) {
	cilk_for(int i=0 ; i<nv ; i++) {
		dist[ii] = MYINFINITY;
	}
}


//--------------------------------------------------------------------------------------------------
// - Validation checking for all the edges to get next vertex frontier

void processedge2(foru *dist, Graph &graph, unsigned iteration, int dst, cilk::reducer_list_append< pair<int,int> > &reduce_list) {

  	foru wt = graph.edgessrcwt[dst];
  	if (wt >= MYINFINITY) return 0;

  	foru dstwt = dist[dst];
  	foru altdist = dist[src] + wt;  

  	if(altdist < dstwt && dst < graph.nnodes){
      	dist[dst] = altdist;
      	reduce_list.push_back(dst);
    }
  
}


//--------------------------------------------------------------------------------------------------
// - Process all the nodes from current vertex frontier
// - Will check all the neighbors of vertices from curr VF 
// - pass those neighbors for validation checking

__device__
unsigned processnode2(foru *dist, Graph &graph, Worklist2 &inwl, Worklist2 &outwl, unsigned iteration){
   	int nn;
 	int nitems;
  	nitems = inwl.size();            // Total number of vertices to be processed

  	cilk::reducer_list_append< pair<int,int> > edgelist;           // To store the neighbor and its src

  	cilk_for(List::const_iterator id = inwl.begin(); id != inwl.end(); id++){
    	int neighborsize = 0;		  // Total no. of neighbors
    	int neighboroffset = 0;		  // Offset of current node in Column array
    	int total_edges = 0;

    	nn = *id;     	
 		neighborsize = graph.getOutDegree(nn);
  		neighboroffset = graph.psrc[graph.srcsrc[nn]];

      	/* if(total_edges) */
      	/* 	  printf("total edges: %d\n", total_edges); */

      	int i;
	  	for(i = 0; i < neighborsize ; i++){
	        edgelist.push_back(make_pair(nn,graph.edgessrcdst[i + neighboroffset]);
	    }

	}

   	List edge_worklist;   // Worklist containing current edge frontier
	edge_worklist = edgelist.get_value();

	Cilk_Reduce reduce_list;

	cilk_for(List::const_iterator id = edge_worklist.begin(); id != edge_worklist.end(); id++){
	   processedge2(dist, graph, iteration, *id, reduce_list)
	}

	outwl = reduce_list.get_value(); // Get next vertex frontier
}



//--------------------------------------------------------------------------------------------------
// - Processing the current vertex frontier
// - If iteration no. is 0 the just push src note to worklist
// - else pass the worklist to the processnode function

void drelax(foru *dist, Graph& graph, List &inwl, List &outwl, int iteration) {
	if(iteration == 0){
		Cilk_Reduce reduce_list;
	    int item = 0;
		reduce_list.push_back(item);
		inwl = reduce_list.get_value();
	}
	else{
	    processnode2(dist, graph, iteration, inwl, outwl);  // Processing nodes from cuurent vertex frontier
	}
}

//--------------------------------------------------------------------------------------------------
// - Display the items present in the worklist

void display_items(List &inwl)
  {
    int count = 0;
    printf("WL: ");
    for(List::const_iterator id = inwl.begin(); id != inwl.end(); id++)
      printf("%d %d, ", count++, *id);

    printf("\n");
  }


//--------------------------------------------------------------------------------------------------
// - Main sssp function
// - Modified Bellman Ford Logic 
// - (No need to ralax all edges in each iteration, 
//    if vertex V has a distance value that has not changed since last time the edges out of V
//    were relaxed, then there is no need to relax the edges out of V second time)

void sssp(foru *dist, Graph &graph)
{
	foru foruzero = 0.0;
	int iteration = 0;     // To keep track of no. of iterations

	double starttime, endtime;
	double runtime;

	
	initialize(dist, graph.nnodes);  // Initializing all the disctances from src to infinity

	dist[0] = 0;  // Initialize src distance to 0

	int nitems = 1;
	
	std::list<int> inwl;   // List containg current vertex frontier
	std::list<int> outwl;  // List containing  next vertex frontier

	//print_array<<<1, graph.nedges + 1>>>((int *) graph.edgessrcdst, graph.nedges + 1);
	//print_texture<<<1, graph.nedges + 1>>>(graph.nedges + 1);
	//return;


	printf("solving.\n");
	printf("starting...\n");
	//printf("worklist size: %d\n", inwl->nitems());
	//printf("WL: 0 0, \n");

	starttime = rtclock();
	drelax(dist, graph, nerr, *inwl, *outwl, 0);

	do {
	    ++iteration;
	
		//printf("ITERATION: %d\n", iteration);
		//display_items(inwl);
		
		drelax(dist, graph, *inwl, *outwl, iteration);
		nitems = outwl.size();

		//remove_dups<<<14, 1024>>>(*outwl, node_owners, gb);
		
		//printf("%d\n", iteration);
		//display_items(outwl);

		//printf("worklist size: %d\n", nitems);
		inwl.clear();
		cilk_for(List::const_iterator id = outwl.begin(); id != outwl.end(); id++){  // Copying the next vertex frontier into current VF
			inwl.push_back(*id);
		}
		outwl.clear();

	} while (nitems > 0);

	endtime = rtclock();

	//printf("\titerations = %d %d.\n", iteration, hwp.iteration);
	//runtime = (1000.0f * (endtime - starttime));
	//printf("\truntime [%s] = %f ms.\n", SSSP_VARIANT, runtime);

	return;
}

