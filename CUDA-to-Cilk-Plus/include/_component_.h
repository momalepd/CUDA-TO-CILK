struct ComponentSpace {
	ComponentSpace(unsigned nelements);
	
	unsigned numberOfElements();
	unsigned numberOfComponents();
	bool isBoss(unsigned element);
	unsigned find(unsigned lelement, bool compresspath = true);
	bool unify(unsigned one, unsigned two);
	void print();
    void copy(ComponentSpace &two);
   
    void dump_to_file(const char *F);
	void allocate();
	void init();

	unsigned nelements;
	unsigned ncomponents,			// number of components.
		 	 *complen, 				// lengths of components.
			 *ele2comp;				// components of elements.
};

//--------------------------------------------------------------------------------------------------------

ComponentSpace::ComponentSpace(unsigned nelements) {
	this->nelements = nelements;

	allocate();
	init();
}

//--------------------------------------------------------------------------------------------------------

void ComponentSpace::dump_to_file(const char *F)
{
  static FILE *f;

  if(!f){
      f = fopen(F, "w");
    }

  int i;
  for(i = 0; i < nelements; i++)
    {
      fprintf(f, "%d %d\n", i, ele2comp[i]);
    }
  fprintf(f, "\n");
}

//--------------------------------------------------------------------------------------------------------

void ComponentSpace::copy(ComponentSpace &two)
{
  two.ncomponents = ncomponents ;
  cilk_for (int id = 0 ; id < nelements ; id++){
  	two.ele2comp[id] = ele2comp[id];
  	two.complen[id] = complen[id];
  }
}

//--------------------------------------------------------------------------------------------------------

void ComponentSpace::print() {
	printf("\t\t-----------------\n");
	for (unsigned id = 0; id < nelements; ++id) {
		printf("\t\t%d -> %d\n", id, ele2comp[id]);
	}	
	printf("\t\t-----------------\n");
}

//--------------------------------------------------------------------------------------------------------

unsigned ComponentSpace::numberOfElements() {
	return nelements;
}

unsigned ComponentSpace::numberOfComponents() {
	return ncomponents;
}

//--------------------------------------------------------------------------------------------------------

void ComponentSpace::allocate() {
	complen = (unsigned *)malloc(nelements * sizeof(unsigned)); 
	ele2comp = (unsigned *)malloc(nelements * sizeof(unsigned)); 
}

//--------------------------------------------------------------------------------------------------------

void ComponentSpace::init() {
	ncomponents = nelements;
	cilk_for (int id = 0 ; id < nelements ; id++) {
		//elements[id] 	= id;
		complen[id]	= 1;
		ele2comp[id] = id;
	}
}

//--------------------------------------------------------------------------------------------------------
// %%

bool ComponentSpace::isBoss(unsigned element) {
  return (ele2comp[element] == element) ? true : false ;
}

unsigned ComponentSpace::find(unsigned lelement, bool compresspath/*= true*/) {
	// do we need to worry about concurrency in this function?
	// for other finds, no synchronization necessary as the data-structure is a tree.
	// for other unifys, synchornization is not required considering that unify is going to affect only bosses, while find is going to affect only non-bosses.
	unsigned element = lelement;
	while (isBoss(element) == false) {
	  element = ele2comp[element];
	}
	if (compresspath) ele2comp[lelement] = element;	// path compression.
	return element;
}

//%%
bool ComponentSpace::unify(unsigned one, unsigned two) {
	// if the client makes sure that one component is going to get unified as a source with another destination only once, then synchronization is unnecessary.
	// while this is true for MST, due to load-balancing in if-block below, a node may be source multiple times.
	// if a component is source in one thread and destination is another, then it is okay for MST.
    do {
      if(!isBoss(one)) return false;
      if(!isBoss(two)) return false;

      unsigned onecomp = one;
      unsigned twocomp = two;
      //unsigned onecomp = find(one, false);
      //unsigned twocomp = find(two, false);

      if (onecomp == twocomp) return false; // "duplicate" edges due to symmetry

		unsigned boss = twocomp;
		unsigned subordinate = onecomp;
		//if (complen[onecomp] > complen[twocomp]) {	// one is larger, make it the representative: can create cycles.
		if (boss < subordinate) {			// break cycles by id.
			boss = onecomp;
			subordinate = twocomp;
		}
		
		// merge subordinate into the boss.
		//ele2comp[subordinate] = boss;

		unsigned oldboss = atomicCAS(&ele2comp[subordinate], subordinate, boss);
		
		if (oldboss != subordinate) {	// someone else updated the boss.
		// we need not restore the ele2comp[subordinate], as union-find ensures correctness and complen of subordinate doesn't matter.
			one = oldboss;
			two = boss;
			return false;
		} else {
			printf("\t\tunifying %d -> %d (%d)\n", subordinate, boss);
			atomicAdd(&complen[boss], complen[subordinate]);
			//complen[boss] += complen[subordinate];
			// complen[subordinate] doesn't matter now, since find() will find its boss.
	
			// a component has reduced.
			unsigned ncomp = atomicSub(ncomponents, 1);
			//atomicDec(ncomponents, nelements);
			printf("\t%d: ncomponents = %d\n", threadIdx.x, ncomp);
			return true;
		}
    
    } while (true);
}
