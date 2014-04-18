/* N-Body simulation Barnes Hut tree based appraoch
	- Implementation of CUDA version of GPU to Cilk Plus on Multicore Processor
	- Lonestargpu benchmark is used to migrate from CUDA to Cilk Plus
	- This is a plain transformation of CUDA BH to Cilk Plus.
*/

#include <cilk/cilk.h>
#include <atomic.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sys/time.h>
#include <cilk/reducer_max.h>
#include <cilk/reducer_min.h>
#include <cilk/reducer_opadd.h>


volatile int step, bottom, maxdepth;
volatile float radius;

//----------------------------------------------------------------------------------------------------
// compute center and radius

void BoundingBoxKernel(int nnodes, int nbodies, int * start, int * child, float * mass, float * posx, float * posy, float * posz)
{
  
  	register float val, maxx_, maxy_, maxz_, minx_, miny_, minz_;
  	cilk::reducer_max<float> maxx, maxy, maxz;
  	cilk::reducer_min<float> minx, miny, minz;


  	cilk_for(int i = 0 ; i < nbodies ; i++) {
  		minx.calc_min(posx[i]);
  		miny.calc_min(posy[i]);
  		minz.calc_min(posz[i]);
  		maxx.calc_max(posx[i]);
  		maxy.calc_max(posy[i]);
  		maxz.calc_max(posz[i]);
  	}

  	maxx_ = maxx.get_value();
  	maxy_ = maxy.get_value();
  	maxz_ = maxz.get_value();
  	minx_ = minx.get_value();
  	miny_ = miny.get_value();
  	minz_ = minz.get_value();

    // compute 'radius'
    val = fmaxf(maxx_ - minx_, maxy_ - miny_);
    radius = fmaxf(val, maxz_ - minz_) * 0.5f;

    // create root node
    k = nnodes;
    bottom = k;

    massd[k] = -1.0f;
    startd[k] = 0;
    posx[k] = (minx_ + maxx_) * 0.5f;
    posy[k] = (miny_ + maxy_) * 0.5f;
    posz[k] = (minz_ + maxz_) * 0.5f;
    k *= 8;
    for (i = 0; i < 8; i++) child[k + i] = -1;

    step++;
}
//------------------------------------------------------------------------------------------------

void Clear1(int nnodes, int nbodies, int * child)  // Store second half of child array to -1
{
  register int top, bot;

  top = 8 * nnodes;
  bot = 8 * nbodies;

  // iterate over all cells from bottom to top
  cilk_for (int k = bot ; k < top ; k++) {
    child[k] = -1;
  }
}


void TreeBuildingKernel(int nnodes, int nbodies, int * child, float * posx, float * posy, float * posz)
{
  	
  	cilk_for (int i = 0 ; i < nbodies ; i++){	

  		int j, depth, localmaxdepth, skip;
  		float x, y, z, r;
  		float px, py, pz;
    	float dx, dy, dz;
    	int ch, n, cell, locked, patch;

		// cache root data
		float radi = radius;
		float rootx = posx[nnodes];
		float rooty = posy[nnodes];
		float rootz = posz[nnodes];

  		int localmaxdepth = 1;
  		int skip = 1;
  		bool change = false;
  		// iterate over all bodies assigned to thread
  		while (false == change) {
    		if (skip != 0) {
      			// new body, so start traversing at root
      			skip = 0;
      			px = posx[i];
      			py = posy[i];
      			pz = posz[i];
      			n = nnodes;
      			depth = 1;
      			r = radi * 0.5f;
      			dx = dy = dz = -r;
      			j = 0;
      			// determine which child to follow
      			if (rootx < px) {j = 1; dx = r;}
      			if (rooty < py) {j |= 2; dy = r;}
      			if (rootz < pz) {j |= 4; dz = r;}
      			x = rootx + dx;
      			y = rooty + dy;
 			    z = rootz + dz;
    		}

    		// follow path to leaf cell
    		ch = child[n*8+j];  // J'th child of root
    		while (ch >= nbodies) {
      			n = ch;
      			depth++;
      			r *= 0.5f;
      			dx = dy = dz = -r;
      			j = 0;
      			// determine which child to follow
      			if (x < px) {j = 1; dx = r;}
      			if (y < py) {j |= 2; dy = r;}
      			if (z < pz) {j |= 4; dz = r;}
      			x += dx;
      			y += dy;
      			z += dz;
      			ch = child[n*8+j];
    		}

    		if (ch != -2) {  // skip if child pointer is locked and try again later
      			locked = n*8+j;
      			if (ch == -1) {
        			if (-1 == atomicCAS((int *)&child[locked], -1, i)) {  // if null, just insert the new body
          				localmaxdepth = max(depth, localmaxdepth);
          				change = true;
        			}
      			} else {  // there already is a body in this position
        			if (ch == atomicCAS((int *)&child[locked], ch, -2)) {  // try to lock
          				patch = -1;
	          			// create new cell(s) and insert the old and new body
	          			do {
	            			depth++;

	            			cell = atomicSub((int *)&bottom, 1) - 1;
				            if (cell <= nbodies) {
	              				//*errd = 1;
	              				bottom = nnodes;
	            			}

	            			if (patch != -1) {
	              				child[n*8+j] = cell;
	            			}
	            			patch = max(patch, cell);

	            			j = 0;
	            			if (x < posx[ch]) j = 1;
	            			if (y < posy[ch]) j |= 2;
	            			if (z < posz[ch]) j |= 4;
	            			child[cell*8+j] = ch;

	            			n = cell;
	            			r *= 0.5f;
	            			dx = dy = dz = -r;
	            			j = 0;
	            			if (x < px) {j = 1; dx = r;}
	            			if (y < py) {j |= 2; dy = r;}
	            			if (z < pz) {j |= 4; dz = r;}
	            			x += dx;
	            			y += dy;
	            			z += dz;

	            			ch = child[n*8+j];
	            		// repeat until the two bodies are different children
	          			} while (ch >= 0);
	          			
	          			child[n*8+j] = i;

	          			localmaxdepth = max(depth, localmaxdepth);
	       			    skip = 2;
	       			    change = true;
        			}
      			}
    		}
    

    		if (skip == 2) {
      			child[locked] = patch;
   			}
  		}
  	// record maximum tree depth
 		atomicMax((int *)&maxdepth, localmaxdepth);
	}
}


void Clear2(int nnodes, int * start, float * mass)
{
  
  // iterate over all cells which are intermediate nodes
 	cilk_for (int k = bottom ; k < nnodes ; k++) {
    	mass[k] = -1.0f;
    	start[k] = -1;
  	}
}

//------------------------------------------------------------------------------------------------
//%% compute center of mass

void SummarizationKernel(const int nnodes, const int nbodies, int * count, const int * child, float * massd, float * posx, float * posy, float * posz)
{
  
  register int i, j, k, ch, inc, cnt, bottom, flag;
  register float m, cm, px, py, pz;
  
  bottom = bottomd;
  inc = blockDim.x * gridDim.x;
  k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
  if (k < bottom) k += inc;

  
  // iterate over all cells assigned to thread
  cilk_for (int k = bottom ; k <= nnodes ; k++) {
  	
  	int i, j, ch, cnt, flag;
  	float m, cm, px, py, pz;
  	int childl[8];
	float massl[8];

    if (mass[k] < 0.0f) {

    	flag = 0;
    	j = 0;

    	while(flag == 0){

    		if (j == 0) {
        		j = 8;
        		for (i = 0; i < 8; i++) {
          			ch = child[k*8+i];
          			childl[i] = ch;  // cache children
          			if ((ch < nbodies) || ((massl[i] = mass[ch]) >= 0.0f)) {
            			j--;
          			}
        		}
      		} else {
        			j = 8;
        			for (i = 0; i < 8; i++) {
          				ch = childl[i];
          				if ((ch < nbodies) || (massl[i] >= 0.0f) || ((massl[i] = mass[ch]) >= 0.0f)) {
            				j--;
          				}
        			}	
      		}

      		if (j == 0) {
        		// all children are ready
		        cm = 0.0f;
		        px = 0.0f;
		        py = 0.0f;
		        pz = 0.0f;
		        cnt = 0;
		        for (i = 0; i < 8; i++) {
		          	ch = childl[i];
		          	if (ch >= 0) {
		            	if (ch >= nbodies) {  // count bodies (needed later)
		              		m = massl[i];
		              		cnt += count[ch];
		            	} else {
		              		m = mass[ch];
		              		cnt++;
		            	}
		            	// add child's contribution
		            	cm += m;
		            	px += posx[ch] * m;
		            	py += posy[ch] * m;
		            	pz += posz[ch] * m;
          			}
        		}
        		count[k] = cnt;
        		m = 1.0f / cm;
		        posx[k] = px * m;
		        posy[k] = py * m;
		        posz[k] = pz * m;
		        flag = 1;
      		}
      	}
    	mass[k] = cm;
    }
  }
}

//------------------------------------------------------------------------------------------------
// Sort Bodies

void SortKernel(int nnodes, int nbodies, int * sort, int * count, int * start, int * child)
{
	// iterate over all cells assigned to thread
	cilk_for (int k = nnodes ; k >= bottom ; k--) {

		int i, j, ch, st,flag;
	    st = start[k];
	    flag = 0;

	    while(flag == 0){

		    if (st >= 0) {
		    	flag = 1;
		      	j = 0;
		      	for (i = 0; i < 8; i++) {
		        	ch = child[k*8+i];
		        	if (ch >= 0) {
		          		if (i != j) {
		            		// move children to front (needed later for speed)
		            		child[k*8+i] = -1;
		            		child[k*8+j] = ch;
		          		}
		         		j++;
		          		if (ch >= nbodies) {
		            		// child is a cell
		            		start[ch] = st;  // set start ID of child
		            		st += count[ch];  // add #bodies in subtree
		          		} else {
		            		// child is a body
		            		sort[st] = ch;  // record body in 'sorted' array
		            		st++;
		          		}	
		        	}
		      	}
		    }
		}
	}
}

//------------------------------------------------------------------------------------------------
// compute force


//------------------------------------------------------------------------------------------------
// advance bodies

void IntegrationKernel(int nbodies, float dtime, float dthf, float * posx, float * posy, float * posz, float *  velxd, float * vely, float * velz, float * accx, float * accy, float * accz) {
	
	// iterate over all bodies assigned to thread
	cilk_for (i = 0; i < nbodies; i ++) {
		float dvelx, dvely, dvelz;
		float velhx, velhy, velhz;
	
	    // integrate
	    dvelx = accx[i] * dthf;
	    dvely = accy[i] * dthf;
	    dvelz = accz[i] * dthf;

	    velhx = velx[i] + dvelx;
	    velhy = vely[i] + dvely;
	    velhz = velz[i] + dvelz;

	    posx[i] += velhx * dtime;
	    posy[i] += velhy * dtime;
	    posz[i] += velhz * dtime;

	    velx[i] = velhx + dvelx;
	    vely[i] = velhy + dvely;
	    velz[i] = velhz + dvelz;
	}
}


//------------------------------------------------------------------------------------------------
// random number generator

#define MULT 1103515245
#define ADD 12345
#define MASK 0x7FFFFFFF
#define TWOTO31 2147483648.0

static int A = 1;
static int B = 0;
static int randx = 1;
static int lastrand;


static void drndset(int seed)
{
   A = 1;
   B = 0;
   randx = (A * seed + B) & MASK;
   A = (MULT * A) & MASK;
   B = (MULT * B + ADD) & MASK;
}


static double drnd()
{
   lastrand = randx;
   randx = (A * randx + B) & MASK;
   return (double)lastrand / TWOTO31;
}


//------------------------------------------------------------------------------------------------

int main(int argc, char *argv[]){

	register int i, run;
	int nnodes, nbodies, step, timesteps;
	register double runtime;
	int error;
	register float dtime, dthf, epssq, itolsq;
	float time, timing[7];
	float *mass, *posx, *posy, *posz, *velx, *vely, *velz;

	int  *sort, *child, *count, *start;
	float *massl;
	float *accx, *accy, *accz;
	float *maxxl, *maxyl, *maxzl;
	float *minxl, *minyl, *minzl;
	register double rsc, vsc, r, v, x, y, z, sq, scale;

	fflush(stdout);
  	if (argc != 3) {
    	fprintf(stderr, "\n");
    	fprintf(stderr, "arguments: number_of_bodies number_of_timesteps\n");
    	exit(-1);
  	}

  	for (run = 0; run < 3; run++) {
	    for (i = 0; i < 7; i++) timing[i] = 0.0f;

	    nbodies = atoi(argv[1]);
	    if (nbodies < 1) {
	  	    fprintf(stderr, "nbodies is too small: %d\n", nbodies);
	    	exit(-1);
	    }
	    if (nbodies > (1 << 30)) {
	      	fprintf(stderr, "nbodies is too large: %d\n", nbodies);
	      	exit(-1);
	    }

	    timesteps = atoi(argv[2]);
	    dtime = 0.025;  dthf = dtime * 0.5f;
    	epssq = 0.05 * 0.05;
    	itolsq = 1.0f / (0.5 * 0.5);

    	// allocate memory

    	if (run == 0) {
		    printf("configuration: %d bodies, %d time steps\n", nbodies, timesteps);

		    mass = (float *)malloc(sizeof(float) * nbodies); if (mass == NULL) {fprintf(stderr, "cannot allocate mass\n");  exit(-1);}
		    posx = (float *)malloc(sizeof(float) * nbodies); if (posx == NULL) {fprintf(stderr, "cannot allocate posx\n");  exit(-1);}
		    posy = (float *)malloc(sizeof(float) * nbodies); if (posy == NULL) {fprintf(stderr, "cannot allocate posy\n");  exit(-1);}
		    posz = (float *)malloc(sizeof(float) * nbodies); if (posz == NULL) {fprintf(stderr, "cannot allocate posz\n");  exit(-1);}
		    velx = (float *)malloc(sizeof(float) * nbodies); if (velx == NULL) {fprintf(stderr, "cannot allocate velx\n");  exit(-1);}
		    vely = (float *)malloc(sizeof(float) * nbodies); if (vely == NULL) {fprintf(stderr, "cannot allocate vely\n");  exit(-1);}
		    velz = (float *)malloc(sizeof(float) * nbodies); if (velz == NULL) {fprintf(stderr, "cannot allocate velz\n");  exit(-1);}

		      
		    child = (int *)malloc(sizeof(int) * (nnodes+1) * 8); if (child == NULL) {fprintf(stderr, "cannot allocate child\n");  exit(-1);}
		    accx = (float *)malloc(sizeof(float) * (nnodes+1));  if (accx == NULL) {fprintf(stderr, "cannot allocate accx\n");  exit(-1);}
		    accy = (float *)malloc(sizeof(float) * (nnodes+1));  if (accy == NULL) {fprintf(stderr, "cannot allocate accy\n");  exit(-1);}
		    accz = (float *)malloc(sizeof(float) * (nnodes+1));  if (accz == NULL) {fprintf(stderr, "cannot allocate accz\n");  exit(-1);}
		    count = (int *)malloc(sizeof(int) * (nnodes+1));     if (count == NULL) {fprintf(stderr, "cannot allocate count\n");  exit(-1);}
		    start = (int *)malloc(sizeof(int) * (nnodes+1));     if (start == NULL) {fprintf(stderr, "cannot allocate start\n");  exit(-1);}
		    sort = (int *)malloc(sizeof(int) * (nnodes+1));      if (sort == NULL) {fprintf(stderr, "cannot allocate sort\n");  exit(-1);}
		    
    	}

    	// generate input

	    drndset(7);
	    rsc = (3 * 3.1415926535897932384626433832795) / 16;
	    vsc = sqrt(1.0 / rsc);
	    cilk_for (i = 0; i < nbodies; i++) {
	      	mass[i] = 1.0 / nbodies;
	      	r = 1.0 / sqrt(pow(drnd()*0.999, -2.0/3.0) - 1);
	      	do {
	        	x = drnd()*2.0 - 1.0;
	        	y = drnd()*2.0 - 1.0;
	        	z = drnd()*2.0 - 1.0;
	        	sq = x*x + y*y + z*z;
	      	} while (sq > 1.0);
	      	scale = rsc * r / sqrt(sq);
	      	posx[i] = x * scale;
	      	posy[i] = y * scale;
	      	posz[i] = z * scale;

	      	do {
	        	x = drnd();
	        	y = drnd() * 0.1;
	      	} while (y > x*x * pow(1 - x*x, 3.5));
	      	v = x * sqrt(2.0 / sqrt(1 + r*r));
	      	do {
	        	x = drnd()*2.0 - 1.0;
	        	y = drnd()*2.0 - 1.0;
	        	z = drnd()*2.0 - 1.0;
	        	sq = x*x + y*y + z*z;
	      	} while (sq > 1.0);
	      	scale = vsc * v / sqrt(sq);
	      	velx[i] = x * scale;
	      	vely[i] = y * scale;
	      	velz[i] = z * scale;
	    }

	    // run timesteps

	    step = -1;
	    maxdepth = 1;

	    for (step = 0; step < timesteps; step++) {
	      
	      	BoundingBoxKernel(nnodes, nbodies, start, child, mass, posx, posy, posz);
	      
	      	Clear1(nnodes, nbodies, child);
	      	TreeBuildingKernel(nnodes, nbodies, child, posx, posy, posz);
	      	Clear2(nnodes, startl, mass);
	    
	 	    SummarizationKernel(nnodes, nbodies, count, child, mass, posx, posy, posz);
	 
	 	    SortKernel(nnodes, nbodies, sort, count, start, child);
	      
	      	ForceCalculationKernel(nnodes, nbodies, dthf, itolsq, epssq, sort, child, mass, posx, posy, posz, velx, vely, velz, accx, accy, accz);
	      
	      cudaEventRecord(start, 0);
	      IntegrationKernel<<<blocks * FACTOR6, THREADS6>>>(nbodies, dtime, dthf, posxl, posyl, poszl, velxl, velyl, velzl, accxl, accyl, acczl);
	      cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
	      timing[6] += time;
	      CudaTest("kernel 6 launch failed");
    	}




  

}

