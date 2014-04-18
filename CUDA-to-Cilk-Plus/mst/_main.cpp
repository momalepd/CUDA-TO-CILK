/* Minimum Spanning Tree
	- Implementation of CUDA version of GPU to Cilk Plus on Multicore Processor
	- Lonestargpu benchmark is used to migrate from CUDA to Cilk Plus
	- This is a plain transformation of CUDA MST to Cilk Plus.
*/

#include "lonestargpu.h"
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>

cilk::reducer_opadd<unsigned int> mstwt (0);
cilk::reducer_opadd<unsigned int> edgecount (0);

//-------------------------------------------------------------------------------------------------------------------------------
void init(Graph graph, ComponentSpace cs, foru *eleminwts, foru *minwtcomponent, unsigned *partners, unsigned *phores, bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
	
	cilk_for(int id = 0 ; id < graph.nnodes ; id++) {
		eleminwts[id] = MYINFINITY;
		minwtcomponent[id] = MYINFINITY;	
		goaheadnodeofcomponent[id] = graph.nnodes;
		phores[id] = 0;
		partners[id] = id;
		processinnextiteration[id] = false;
	}
}

//-------------------------------------------------------------------------------------------------------------------------------

void findelemin(Graph graph, ComponentSpace cs, foru *eleminwts, foru *minwtcomponent, unsigned *partners, unsigned *phore, bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
	
	cilk_for (int = id ; id < graph.nnodes ; id++) {
		// if I have a cross-component edge,
		// find my minimum wt cross-component edge,
		// inform my boss about this edge e (atomicMin).
		unsigned src = id;
		unsigned srcboss = cs.find(src);
		unsigned dstboss = graph.nnodes;
		foru minwt = MYINFINITY;
		unsigned degree = graph.getOutDegree(src);

		for (unsigned ii = 0; ii < degree; ++ii) {  // Checking neighbours of each indivisual neighbours
			foru wt = graph.getWeight(src, ii);
			if (wt < minwt) {
				unsigned dst = graph.getDestination(src, ii);
				unsigned tempdstboss = cs.find(dst);
				if (srcboss != tempdstboss) {	// cross-component edge.
					minwt = wt;
					dstboss = tempdstboss;
				}
			}
		}

		printf("\tminwt[%d] = %d\n", id, minwt);
		eleminwts[id] = minwt;
		partners[id] = dstboss;

		if (minwt < minwtcomponent[srcboss] && srcboss != dstboss) {
			//? inform boss.
			minwtcomponent[srcboss], minwt;
		}
	}
}

//-------------------------------------------------------------------------------------------------------------------------------

void findelemin2(Graph graph, ComponentSpace cs, foru *eleminwts, foru *minwtcomponent, unsigned *partners, unsigned *phore, bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
	
	cilk_for (int id = 0 ; id < graph.nnodes ; id++) {
		unsigned src = id;
		unsigned srcboss = cs.find(src);

		if(eleminwts[id] == minwtcomponent[srcboss] && srcboss != partners[id] && partners[id] != graph.nnodes){
		    unsigned degree = graph.getOutDegree(src);
		    for (unsigned ii = 0; ii < degree; ++ii) {
		      	foru wt = graph.getWeight(src, ii);
		      	if (wt == eleminwts[id]) {
					unsigned dst = graph.getDestination(src, ii);
					unsigned tempdstboss = cs.find(dst);
					if (tempdstboss == partners[id]) {	// cross-component edge.
			  			//atomicMin(&goaheadnodeofcomponent[srcboss], id);
			  
			  			if(atomicCAS(&goaheadnodeofcomponent[srcboss], graph.nnodes, id) == graph.nnodes){
			 	     		//printf("%d: adding %d\n", id, eleminwts[id]);
			      			//atomicAdd(wt2, eleminwts[id]);
			    		}
					}	
		      	}
		    }
		}
	}
}


//-------------------------------------------------------------------------------------------------------------------------------

void verify_min_elem(Graph graph, ComponentSpace cs, foru *eleminwts, foru *minwtcomponent, unsigned *partners, unsigned *phore, bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
	
	cilk_for (int id = 0 ; id < graph.nnodes ; id++) {
	  	if(cs.isBoss(id)){
	    	if(goaheadnodeofcomponent[id] != graph.nnodes){
	     		unsigned minwt_node = goaheadnodeofcomponent[id];
			    unsigned degree = graph.getOutDegree(minwt_node);
	      		foru minwt = minwtcomponent[id];

	      		if(minwt != MYINFINITY){
					bool minwt_found = false;
	      			//printf("%d: looking at %d def %d minwt %d\n", id, minwt_node, degree, minwt);
	      			for (unsigned ii = 0; ii < degree; ++ii) {
						foru wt = graph.getWeight(minwt_node, ii);
						//printf("%d: looking at %d edge %d wt %d (%d)\n", id, minwt_node, ii, wt, minwt);

						if (wt == minwt) {
		  					minwt_found = true;
		  					unsigned dst = graph.getDestination(minwt_node, ii);
		  					unsigned tempdstboss = cs.find(dst);
		  					if(tempdstboss == partners[minwt_node] && tempdstboss != id){
		      					processinnextiteration[minwt_node] = true;
		      					//printf("%d okay!\n", id);
		    				}
						}
	      			}	
	      		}
	    	}
		}
	}
}

//-------------------------------------------------------------------------------------------------------------------------------

void elim_dups(Graph graph, ComponentSpace cs, foru *eleminwts, foru *minwtcomponent, unsigned *partners, unsigned *phore, bool *processinnextiteration, unsigned *goaheadnodeofcomponent) {
	
	cilk_for (int id = 0 ; id < graph.nnodes ; id++) {
	  	if(processinnextiteration[id]){
	      	unsigned srcc = cs.find(id);
	      	unsigned dstc = partners[id];
	      
	      	if(minwtcomponent[dstc] == eleminwts[id]){
		  		if(id < goaheadnodeofcomponent[dstc]){
		      		processinnextiteration[id] = false;
		      		//printf("duplicate!\n");
		    	}
			}
	    }
	}
}
//-------------------------------------------------------------------------------------------------------------------------------

void dfindcompmintwo(Graph graph, ComponentSpace csw, foru *eleminwts, foru *minwtcomponent, unsigned *partners, unsigned *phores, bool *processinnextiteration, unsigned *goaheadnodeofcomponent, bool &repeat) {
	unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned id, nthreads = blockDim.x * gridDim.x;
	if (inpid < graph.nnodes) id = inpid;

	unsigned up = (graph.nnodes + nthreads - 1) / nthreads * nthreads;
	unsigned srcboss, dstboss;


	cilk_for(int id = 0; id < graph.nnodes; id ++) {
	 	if(id < graph.nnodes && processinnextiteration[id]){
	      	srcboss = csw.find(id);
	      	dstboss = csw.find(partners[id]);
	    }
	 
	 	if (id < graph.nnodes && processinnextiteration[id] && srcboss != dstboss) {
	    	//printf("trying unify id=%d (%d -> %d)\n", id, srcboss, dstboss);

	    	if (csw.unify(srcboss, dstboss)) {
	      		mstwt += eleminwts[id];
	      		edgecount += 1;
	      		//printf("u %d -> %d (%d)\n", srcboss, dstboss, eleminwts[id]);
	      		processinnextiteration[id] = false;
	      		eleminwts[id] = MYINFINITY;	// mark end of processing to avoid getting repeated.
	    	}
	    	else {
	      		repeat = true;
	    	}
	    	printf("\tcomp[%d] = %d.\n", srcboss, csw.find(srcboss));
	  	}

	}
}

//-------------------------------------------------------------------------------------------------------------------------------


int main(int argc, char *argv[]) {
  int iteration = 0;
  Graph graph;

  unsigned *partners, *phores;
  foru *eleminwts, *minwtcomponent;
  bool *processinnextiteration;
  unsigned *goaheadnodeofcomponent;

  double starttime, endtime;

  if (argc != 2) {
    printf("Usage: %s <graph>\n", argv[0]);
    exit(1);
  }

  graph.read(argv[1]);

  ComponentSpace cs(graph.nnodes);

  // Allocating memory
  eleminwts = (foru *)malloc(graph.nnodes * sizeof(foru));
  minwtcomponent = (foru *)malloc(graph.nnodes * sizeof(foru));
  partners = (unsigned *)malloc(graph.nnodes * sizeof(unsigned));
  phores = (unsigned *)malloc(graph.nnodes * sizeof(unsigned));
  processinnextiteration = (bool *)malloc(graph.nnodes * sizeof(bool));
  goaheadnodeofcomponent = (unsigned *)malloc(graph.nnodes * sizeof(unsigned));


  unsigned prevncomponents, currncomponents = graph.nnodes;
  bool repeat = false;

  printf("finding mst.\n");
  
  starttime = rtclock();
  do {
    	++iteration;
    	prevncomponents = currncomponents;
    	init(mstwt, graph, cs, eleminwts, minwtcomponent, partners, phores, processinnextiteration, goaheadnodeofcomponent);
    	//printf("0 %d\n", cs.numberOfComponentsHost());
    	findelemin(mstwt, graph, cs, eleminwts, minwtcomponent, partners, phores, processinnextiteration, goaheadnodeofcomponent);
    	findelemin2(mstwt, graph, cs, eleminwts, minwtcomponent, partners, phores, processinnextiteration, goaheadnodeofcomponent);
    	verify_min_elem(mstwt, graph, cs, eleminwts, minwtcomponent, partners, phores, processinnextiteration, goaheadnodeofcomponent);
    	//elim_dups(mstwt, graph, cs, eleminwts, minwtcomponent, partners, phores, processinnextiteration, goaheadnodeofcomponent);

    	do {
      		repeat = false;
      		dfindcompmintwo <<<nSM * compmintwo_res, 384>>> (mstwt, graph, cs, eleminwts, minwtcomponent, partners, phores, processinnextiteration, goaheadnodeofcomponent,repeat, gedgecount);
	    } while (repeat); // only required for quicker convergence

	    currncomponents = cs.numberOfComponentsHost();
    	printf("\titeration %d, number of components = %d (%d), mstwt = %u mstedges = %u\n", iteration, currncomponents, prevncomponents, mstwt.get_value(), edgecount.get_value());
  } while (currncomponents != prevncomponents);
  endtime = rtclock();
	
  printf("\tmstwt = %u, iterations = %d.\n", mstwt.get_value(), iteration);
  printf("\t%s result: weight: %u, components: %u, edges: %u\n", argv[1], mstwt.get_value(), currncomponents, edgecount.get_value());
  printf("\truntime [mst] = %f ms.\n", 1000 * (endtime - starttime));

  // cleanup left to the OS.

  return 0;
}
