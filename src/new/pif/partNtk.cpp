
#include "partNtk.h"
namespace ymc {

PartNtk::~PartNtk()
{
	for (auto pNtk : m_vSubNtks)
	{
        Abc_NtkCleanMarkA(pNtk);
		Abc_NtkDelete(pNtk);
	}
	for (auto pNtk : m_vSubNtksMapped)
	{
        Abc_NtkCleanMarkA(pNtk);
		Abc_NtkDelete(pNtk);
	}
	/*
	for (auto pIfPars : m_vpIfPars)
	{
		assert(pIfPars);
		If_DsdManFree(pIfPars->pDsdMan, 0);
		delete pIfPars;
	}
	*/
}

void PartNtk::setIfPars(If_Par_t* pIfPars, int threadId = 0)
{
	memset(pIfPars, 0, sizeof(If_Par_t));
	pIfPars->bIsPif = 1;
	pIfPars->pDsdMan = If_DsdManLoad(m_pDsdLibFile);
	pIfPars->piMaxReqTime = &m_iMaxReqTime;
	pIfPars->factor = 1.0;
	pIfPars->iThreadId = threadId;
    If_DsdManSetNewAsUseless(pIfPars->pDsdMan);
	pIfPars->nLutSize = If_DsdManVarNum(pIfPars->pDsdMan);
	pIfPars->nCutsMax = 8; 
	pIfPars->nFlowIters = 1; 
	pIfPars->nAreaIters = 2; 
	pIfPars->DelayTarget = -1;
	pIfPars->Epsilon = static_cast<float>(0.005);
	pIfPars->fPreprocess = 1;
	pIfPars->fEdge = 1;
	pIfPars->fCutMin = 1;
	pIfPars->fUseDsd = 1;
	pIfPars->fUseDsdTune = 1;
	pIfPars->fTruth = 1;
	return;
}

void PartNtk::init()
{
	if(m_nParts > 0)
	{
		m_vSubNtks.resize(m_nParts);
		m_vSubNtksMapped.resize(m_nParts);
	}
}

void PartNtk::startThread()
{
	yassert(m_nParts == m_vSubNtks.size());
	struct timeval t1,t2;
	double time;
	gettimeofday(&t1, NULL);

	vector<thread> vThreads;
	vThreads.reserve(m_nParts);
	for (int i = 0; i < m_nParts; i++)
		vThreads.push_back(thread(threadWrapper, static_cast<void*>(this), i));
	for (auto iter = vThreads.begin(); iter != vThreads.end(); ++iter)
		iter->join();

	gettimeofday(&t2, NULL);
   	time = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
   	printf("startThread spent time: %f\n", time);
}

void PartNtk::threadWrapper(void* pThis, int id)
{
	PartNtk* pPN = static_cast<PartNtk*>(pThis);
	pPN->threadKernel(id);
}

void PartNtk::threadKernel(int id)
{
	ylog("Now in thread %d with %d objs, %d PIs, %d POs\n", id, m_vSubNtks[id]->nObjs, m_vSubNtks[id]->vPis->nSize, m_vSubNtks[id]->vPos->nSize);
	struct timeval t1,t2;
	double time;
	gettimeofday(&t1, NULL);

	If_Par_t ifPars;
	If_Par_t* pIfPars = &ifPars;
	setIfPars(pIfPars, id);
	m_vSubNtksMapped[id] = Abc_NtkIf(m_vSubNtks[id], pIfPars);
	If_DsdManFree(pIfPars->pDsdMan, 0);

	gettimeofday(&t2, NULL);
   	time = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
   	printf("Thread %d spent time: %f\n", id, time);
}

void PartNtk::serialMap()
{
	struct timeval t1,t2;
	double time;
	ylog("Now in serialMap()\n");
	for(int i = 0; i < m_nParts; i++)
	{
		gettimeofday(&t1, NULL);

		If_Par_t ifPars;
		If_Par_t* pIfPars = &ifPars;
		setIfPars(pIfPars);
		m_vSubNtksMapped[i] = Abc_NtkIf(m_vSubNtks[i], pIfPars);
		If_DsdManFree(pIfPars->pDsdMan, 0);

		gettimeofday(&t2, NULL);
	   	time = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
	   	printf("thread %d spent time: %f\n", i, time);
	}
}

void PartNtk::partOriginNtk()
{
	struct timeval t1,t2;
	double time;
	gettimeofday(&t1, NULL);

	if(m_nParts == 0)
	{
        	ylog("No nParts specified. Adaptive partitioning routine is used\n");
		MetisGraph graph(m_pOriginNtk, 0);
		graph.set_sCluster(m_sCluster);
		MetisAig aig;
		aig.bindGraph(&graph);
		//aig.check();
		aig.parseAig();
		m_nParts = aig.partitionAig(); //The Graph is partitioned!
		ylog("After adaptive routine, nParts = %d\n", m_nParts);
		init(); //allocate memory for member vectors
		graph.createSubNtksFromPartition(m_vSubNtks);
	}
	/* comment metis by zli
	else
	{
		ylog("nParts is set to %d\n", m_nParts);
		MetisGraph graph(m_pOriginNtk, 1);
		graph.partGraphByMetis(m_nParts);
		graph.createSubNtksFromPartition(m_vSubNtks);
	}
	*/


	gettimeofday(&t2, NULL);
   	time = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
   	printf("partOriginNtk spent time: %f\n", time);

	//Output the partitioned network
	Abc_NtkWriteVerilog();
}


void PartNtk::mergeMappedSubNtk()
{
	yassert(m_nParts == m_vSubNtks.size());
	yassert(m_vSubNtks.size() == m_vSubNtksMapped.size());
	//Abc_NtkWriteMappedBlif();
	struct timeval t1,t2;
	double time;
	gettimeofday(&t1, NULL);

	Vec_Ptr_t* pSubNtks;
	Vec_Ptr_t* pSubNtksMapped;
    pSubNtks = Vec_PtrAlloc( m_vSubNtks.size() );
    pSubNtksMapped = Vec_PtrAlloc( m_vSubNtksMapped.size() );
	for(auto pNtk : m_vSubNtks)
        Vec_PtrPush( pSubNtks, pNtk);
	for(auto pNtk : m_vSubNtksMapped)
        Vec_PtrPush( pSubNtksMapped, pNtk);

	m_pMappedNtk = Abc_NtkMerge(m_pOriginNtk, pSubNtks, pSubNtksMapped);
	//m_pMappedNtk = Abc_NtkMerge(m_pOriginNtk, pSubNtks, pSubNtks);


	Vec_PtrFree(pSubNtks);
	Vec_PtrFree(pSubNtksMapped);

	gettimeofday(&t2, NULL);
   	time = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
   	printf("mergeMappedSubNtk spent time: %f\n", time);
}

void PartNtk::debug()
{
	printf("Now in debug()\n");
	If_Par_t* pIfPars = new If_Par_t;
	setIfPars(pIfPars);
	//pIfPars->bIsPif = 0;
	Abc_Ntk_t* pNtkMapped = Abc_NtkIf(m_pOriginNtk, pIfPars);
	printf("single if\n");
	Abc_NtkPrintStats(pNtkMapped, 0,0,0,0,0,0,0,0,0,0); 
	/*
	Abc_Obj_t* pCoInNtkMapped = Abc_NtkPickCriticalPo(pNtkMapped);
	if(pCoInNtkMapped == NULL) return;
	int iObjId = Nm_ManFindIdByName(m_pOriginNtk->pManName, Abc_ObjName(pCoInNtkMapped), Abc_ObjType(pCoInNtkMapped) );
	Abc_Obj_t* pCoInNtkOrigin = Abc_NtkObj(m_pOriginNtk, iObjId);
	assert(Abc_ObjIsCo(pCoInNtkOrigin));
	Abc_Ntk_t* pNtkCritical = Abc_NtkExtractCriticalPath(m_pOriginNtk, pCoInNtkOrigin);
	Abc_Ntk_t* pNtkCriticalMapped = Abc_NtkIf(pNtkCritical, m_vpIfPars[0]);
	Abc_NtkPrintStats(pNtkCriticalMapped, 0,0,0,0,0,0,0,0,0,0); 
	*/
}


void PartNtk::Abc_NtkWriteBlif(){
	
	std::cout<<"Now in write blif!"<<endl;
	const char* outputFolder = m_dirName;
	
	int i = 0;
	//char *pFileNameOut = "output.blif";
	for(auto pNtk : m_vSubNtks){

		char filename[256];
		snprintf(filename, sizeof(filename),"%s/network_%d.blif",outputFolder, i);
		printf("Writing file: %s\n", filename);

		Abc_Ntk_t* pNtkNew = Abc_NtkToNetlist(pNtk);
		//cout<<"Now in iteration!"<<endl;

		Io_WriteBlif(pNtkNew,filename,1,0,0);
		i++;
		
	}
}

void PartNtk::Abc_NtkWriteAIG(){
	std::cout<<"Now in write AIG!"<<endl;
	const char* outputFolder = m_dirName;
	
	int i = 0;
	for(auto pNtk : m_vSubNtks){

		char filename[256];
		snprintf(filename, sizeof(filename),"%s/network_%d.aig",outputFolder, i);
		printf("Writing file: %s\n", filename);

		Abc_Ntk_t* pNtkNew = Abc_NtkToNetlist(pNtk);
		Abc_NtkSetName(pNtkNew, const_cast<char*>(("network_" + to_string(i)).c_str()) );
		//cout<<"Now in iteration!"<<endl;
        Io_WriteAiger( pNtk, filename, 1, 0, 0 );
		i++;	
	}
}

void PartNtk::Abc_NtkWriteVerilog(){
	
	std::cout<<"Now in write verilog!"<<endl;
	const char* outputFolder = m_dirName;
	
	int i = 0;
	for(auto pNtk : m_vSubNtks){

		char filename[256];
		snprintf(filename, sizeof(filename),"%s/network_%d.v",outputFolder, i);
		printf("Writing file: %s\n", filename);

		Abc_Ntk_t* pNtkNew = Abc_NtkToNetlist(pNtk);
		Abc_NtkSetName(pNtkNew, const_cast<char*>(("network_" + to_string(i)).c_str()) );
		//cout<<"Now in iteration!"<<endl;

		if ( !Abc_NtkHasAig(pNtkNew) && !Abc_NtkHasMapping(pNtkNew) )
            Abc_NtkToAig( pNtkNew );
        Io_WriteVerilog( pNtkNew, filename, 0 );
		i++;
		
	}
}

void PartNtk::Abc_NtkWriteMappedBlif(){
	std::cout<<"Now in write mapped blif!"<<endl;
	const char* outputFolder = m_dirName;

	int i = 0;
	//char *pFileNameOut = "output.blif";
	for(auto pNtk : m_vSubNtksMapped){

		char filename[256];
		snprintf(filename, sizeof(filename),"%s/Mappednetwork_%d.blif",outputFolder, i);
		printf("Writing file: %s\n", filename);

		Abc_Ntk_t* pNtkNew = Abc_NtkToNetlist(pNtk);
		//cout<<"Now in iteration!"<<endl;

		Io_WriteBlif(pNtkNew,filename,1,0,0);
		i++;
		
	}
}






} //for namespace