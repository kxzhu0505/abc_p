
#include "c_wrapper.h"

namespace ymc {

int hello()
{
	printf("Now in ymc_hello()!\n");
	return 0;
}

int try_metis()
{
	std::cout << "Now in ymc_try_metis()! (printed by iostream)" << std::endl;	
	vector<idx_t> vadj = {0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44};
	vector<idx_t> vadjncy = {1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13,0};
	vector<idx_t> vPart(16);
	idx_t edgecut = 0;
	idx_t nvtx = 15;
	idx_t ncon = 1;
	idx_t nparts = 3;

	METIS_PartGraphKway(&nvtx, &ncon, &vadj[0], &vadjncy[0], NULL, NULL, NULL, &nparts, NULL, NULL, NULL, &edgecut, &vPart[0]);
	cout << "Metis result: " << endl;
	for(int i = 0; i < 15; i++)
		cout << vPart[i] << ' ';
	cout << "\nedgecut = " << edgecut << endl;

	return 0;
}

int test_yaig()
{
	printf("Now in test_yaig()!\n");
    int32_t node[100];

    MetisAig y;
    node[0] = y.addPi();
    node[1] = y.addPi();
    node[2] = y.addPi();
    node[3] = y.addAndNode(node[0], 0, node[1], 1);
    node[4] = y.addAndNode(node[1], 1, node[2], 0);
    node[5] = y.addAndNode(node[3], 0, node[4], 1);
    node[6] = y.addPo(node[5], 0);
	node[7] = y.addPi();
	node[8] = y.addAndNode(node[4], 0, node[7], 1);
	node[9] = y.addPo(node[8], 0);
	node[10] = y.addPi();
	node[11] = y.addPi();
	node[12] = y.addPi();
	node[13] = y.addAndNode(node[10], 0, node[11], 0);
	node[14] = y.addAndNode(node[13], 0, node[12], 0);
	node[15] = y.addPo(node[14], 0);
	y.computeAllLevel();
	y.parseAig();

/*
    y.check();
    shared_ptr<Graph> spg = y.toGraph();
    spg->print();
    spg->check();

    Yaig y2;
    y2.addSubGraph(*spg);
    y2.check();
    spg = y2.toGraph();
    spg->print();
    spg->check();
*/


	return 0;
}

Abc_Ntk_t* pif(Abc_Ntk_t* pNtk, uint32_t nParts, char* libFileName)
{
#ifdef PIF_MULTITHREAD
	ylog("PIF_MULTITHREAD is on\n");
#else
	ylog("PIF_MULTITHREAD is off\n");
#endif
	struct timeval t1,t2;
	double time;
	gettimeofday(&t1, NULL);

	shared_ptr<PartNtk> spPN = make_shared<PartNtk>(pNtk, nParts, libFileName);

	spPN->partOriginNtk();
	spPN->startThread();
//	spPN->serialMap();
	spPN->mergeMappedSubNtk();

	gettimeofday(&t2, NULL);
   	time = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
   	printf("pif spent time: %f\n", time);

	return spPN->getResNtk();
}

} //for namespace

extern "C" int ymc_hello_wrapper()
{
	return ymc::hello();
}

extern "C" int ymc_try_metis_wrapper()
{
	return ymc::try_metis();
}

extern "C" int ymc_test_yaig_wrapper()
{
	return ymc::test_yaig();
}

extern "C" Abc_Ntk_t* ymc_pif_wrapper(Abc_Ntk_t* pNtk, uint32_t nParts, char* libFileName)
{
	return ymc::pif(pNtk, nParts, libFileName);
}
