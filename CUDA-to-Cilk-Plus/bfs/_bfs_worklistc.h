#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <list>


#define cilk::reducer_list_append<int> Cilk_Reduce 
#define std::list<int> List


//---------------------------------------------------------------------------------------
// - Function to initialize all distances to INFINITY

void initialize(foru *dist, unsigned int nv) {
	
	cilk_for(int i=0 ; i<nv ; i++) {
		dist[ii] = MYINFINITY;
	}
}

//---------------------------------------------------------------------------------------
// - Validation checking for all the edges to get next vertex frontier

void processedge2(foru *dist, Graph &graph, unsigned iteration, int dst, Cilk_Reduce &reduce_list) {

  if(dist[dst] == MYINFINITY && dst < graph.nnodes){
      dist[dst] = iteration;
      reduce_list.push_back(dst);
    }
  
}

//---------------------------------------------------------------------------------------
// - Processing all the vertices present in current vertex frontier

unsigned processnode2(foru *dist, Graph &graph, unsigned iteration, List &inwl, List &outwl) {

   	int nn;
   	int nitems;
   	nitems = inwl.size();    // Total number of vertices to be processed

   	/*std::vector<int> worklist;
   	cilk_for(std::list<char>::const_iterator i = values.begin(); i != values.end(); i++)
    {
        worklist.push_back(*i);
    }*/
  
  	Cilk_Reduce edgelist;  // To store the edges of current vertex frontier

   	cilk_for(List::const_iterator id = inwl.begin(); id != inwl.end(); id++){

      	int neighborsize = 0;
      	int neighboroffset = 0;
	    int total_edges = 0;

 		nn = *id;     	
 		neighborsize = graph.getOutDegree(nn);
  		neighboroffset = graph.psrc[graph.srcsrc[nn]];
	    	
	  	int i;
	  	for(i = 0; i < neighborsize ; i++){
	      edgelist.push_back() = graph.edgessrcdst[i + neighboroffset];
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


//---------------------------------------------------------------------------------------


void drelax(foru *dist, Graph& graph, int iteration, List &inwl, List &outwl) {

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


//---------------------------------------------------------------------------------------
// - main bfs functions


void bfs(Graph &graph, foru *dist){}

	foru foruzero = 0;  // foru defined in common.h
	int iteration = 0;  // To keep track of bfs iterations and also to store distance value from src


	initialize(dist, graph.nnodes); // Initialize all the distances to infinity

	int nitems = 1;

	std::list<int> inwl;   // List containg current vertex frontier
	std::list<int> outwl;  // List containing  next vertex frontier

	printf("solving.\n");
	printf("starting...\n");

	drelax(dist, graph,0,inwl,outwl);  // Push src node to worklist
	values = inwl.get_value();
	

	while(nitems > 0){

		++iteration;
		drelax(dist, graph, iteration, inwl, outwl); // Start processing the current vertex frontier

		inwl.clear();
		cilk_for(List::const_iterator id = outwl.begin(); id != outwl.end(); id++){
			inwl.push_back(*id);
		}

		outwl.clear();
		nitems = inwl.size();
	}

	printf("Finished...\n");


}