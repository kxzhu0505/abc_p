#pragma once

#include "yaig.h"


extern "C" {
#include "base/abc/abc.h"
#include "map/if/if.h"
#include <sys/time.h>
extern Abc_Ntk_t * Abc_NtkIf(Abc_Ntk_t * pNtk, If_Par_t * pPars);
}
#include "new/pif/partNtkFuncs.h"

using std::vector;

namespace ymc {

//class for parallel mapping implemented by graph partitioning.
class PartNtk
{
public:
	~PartNtk();
	PartNtk(Abc_Ntk_t* pNtkOrigin, uint32_t nParts, char* libfile): 
		m_nParts(nParts), m_pOriginNtk(pNtkOrigin){m_iMaxReqTime = 0; strcpy(m_pDsdLibFile, libfile); init();};
	void init();
	void setIfPars(If_Par_t* pIfPars, int threadId);
	void setOriginNtk(Abc_Ntk_t* pNtk){ m_pOriginNtk = pNtk; m_pMappedNtk = NULL; }
	Abc_Ntk_t* getResNtk(){ return m_pMappedNtk; }
	uint32_t getNParts(){ return m_nParts; }
	char* getLibFileName() { return m_pDsdLibFile; }

	void partOriginNtk(); //generate m_vSubNtks
	void mergeMappedSubNtk();
	void serialMap();

	static void threadWrapper(void*, int id);
	void threadKernel(int id);
	void startThread();

	void debug();

private:
	uint32_t m_nParts;
	Abc_Ntk_t* m_pOriginNtk;
	Abc_Ntk_t* m_pMappedNtk; 
	int m_iMaxReqTime;
	//vector<If_Par_t*> m_vpIfPars;
	vector<Abc_Ntk_t*> m_vSubNtks; //remember to dealloc
	vector<Abc_Ntk_t*> m_vSubNtksMapped; //remember to dealloc
	char m_pDsdLibFile[100];
};
















} //for namespace
