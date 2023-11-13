#include "partNtkFuncs.h"

namespace ymc{

MetisGraph::MetisGraph(Abc_Ntk_t* pNtk, bool fMetis)
{
	yassert(Abc_NtkIsStrash(pNtk) && Abc_NtkHasAig(pNtk));
	int32_t nNodes = Abc_NtkObjNum(pNtk) - Abc_NtkBoxNum(pNtk) - 1;

	metisGraphInit(nNodes, fMetis);
	createGraphFromNtk(pNtk, fMetis);
}

void MetisGraph::metisGraphInit(int32_t nNodes, bool fMetis)
{
	m_nNodes = nNodes;
	m_nCon = 1;
	m_nPis = 0;
	m_vNodeIndices.reserve(nNodes + 2);
	m_vNodeIsPi.reserve(nNodes + 1);
	m_nParts = 0;
	m_nCutEdges = 0;
	m_vPartition.resize(nNodes + 1, -1);
	m_vpObjs.reserve(nNodes + 1);
	m_vpObjsNew.reserve(nNodes + 1);
	m_vNodeIndices.push_back(0);

	if(fMetis)
	{
		m_vNodeWeights.reserve(nNodes + 1);
		m_vEdgeWeights.reserve(4*nNodes);
	}
	
	m_vEdges.reserve(4*nNodes);
	m_vEdgeIsC.reserve(4*nNodes);
}

int MetisGraph::computeNodeWeightByLevel(int leftChildLevel, int rightChildLevel)
{
	int res = (leftChildLevel>rightChildLevel)?leftChildLevel:rightChildLevel;
	if(res > MAX_NODE_WEIGHT - 1)
		res = MAX_NODE_WEIGHT - 1;
	//return 100;
	return res + 1;
}

int MetisGraph::computeNodeWeightByDegree(int degree)
{
	int res = degree;
	if(res > MAX_NODE_WEIGHT)
		res = MAX_NODE_WEIGHT;
	return res;
}

int MetisGraph::computeEdgeWeight(int faninLevel, int faninReqTime, int level, int reqTime)
{
	if(faninLevel == 0)
		return 1;
	if((reqTime - level < WALL_THRESHOLD) || (faninReqTime - faninLevel < WALL_THRESHOLD))
		return WALL_EDGE_WEIGHT;
	int res = MAX_EDGE_WEIGHT_FOR_NODE;
	int slack = level - 1 - faninLevel;
	if(res - slack*10 > MIN_EDGE_WEIGHT_FOR_NODE)
		res -= slack*10;
	else
		res = MIN_EDGE_WEIGHT_FOR_NODE;
	return res;
}

void MetisGraph::setEdgeWeight(int fanin, int fanout, int weight)
{
	yassert(fanin < fanout);
	yassert(fanout < m_nNodes);
	yassert(getNodeFanin0Id(fanout) == fanin || getNodeFanin1Id(fanout) == fanin);
	int nEdges = getNodeEdgeNum(fanin);
	int iEdge = m_vNodeIndices[fanin];
	bool status = 0;
	for(int i = 0; i < nEdges; i++)
	{
		if(m_vEdges[iEdge + i] == fanout)
		{
			m_vEdgeWeights[iEdge + i] = weight;
			status = 1;
			break;
		}
	}
	yassert(status);
	iEdge = m_vNodeIndices[fanout];
	if(getNodeFanin0Id(fanout) == fanin)
		m_vEdgeWeights[iEdge] = weight;
	else if(getNodeFanin1Id(fanout) == fanin)
		m_vEdgeWeights[iEdge + 1] = weight;
	else
		yassert(0);
}

int MetisGraph::createGraphFromNtk(Abc_Ntk_t* pNtk, bool fMetis)
{
	Abc_Obj_t *pObj;
	Abc_Obj_t *pFanin0, *pFanin1, *pFanout;
	Vec_Ptr_t *vNodes;
	int iNode = 0, iEdgeBeg = 0, iEdgeEnd = 0;
	int iEdge;
	int i, j, weight;
	int maxLevel;
	Abc_NtkFillTemp(pNtk);
	map<int, unsigned int> nodeId2ReqTime;

	m_pNtkOrigin = pNtk;

	if(fMetis)
	{
		maxLevel = Abc_NtkLevelReverse(pNtk);
		Abc_NtkForEachNode(pNtk, pObj, i)
			nodeId2ReqTime.emplace(make_pair(pObj->Id, maxLevel + 1 - pObj->Level)); //BUG?: 2021/6/12: maxLevel - level!!!
		Abc_NtkLevel(pNtk); //mark level for each obj in pNtk
	}



	//Const1
	pObj = Abc_AigConst1(pNtk);
	iEdgeEnd += Abc_ObjFanoutNum(pObj);
	if (iEdgeEnd > iEdgeBeg) //check if the constant node is isolated
	{
		//printf("Const1 node added to the network.\n");
		pObj->iTemp = iNode++; //the corresponding index of the node in Metis Graph
		m_vNodeIsPi.push_back(true);
		m_vpObjs.push_back(pObj);
		m_vpObjsNew.push_back(NULL);
		m_vNodeIndices.push_back(iEdgeEnd);
		iEdgeBeg = iEdgeEnd;
		m_nNodes++;
		m_nPis++;

		if(fMetis)
			m_vNodeWeights.push_back(1);
	}
	//CIs
	Abc_NtkForEachCi(pNtk, pObj, i)
	{
		iEdgeEnd += Abc_ObjFanoutNum(pObj);
		if (iEdgeEnd > iEdgeBeg)
		{
			pObj->iTemp = iNode++;
			m_vNodeIsPi.push_back(true);
			m_vpObjs.push_back(pObj);
			m_vpObjsNew.push_back(NULL);
			m_vNodeIndices.push_back(iEdgeEnd);
			iEdgeBeg = iEdgeEnd;
			m_nPis++;

			if(fMetis)
				m_vNodeWeights.push_back(1);
		}
		else
		{
			//printf("Dangling PI detected.\n");
			m_nNodes--;
		}
	}
	//internal nodes
	vNodes = Abc_NtkDfsIter(pNtk, 0);
	Vec_PtrForEachEntry(Abc_Obj_t *, vNodes, pObj, i)
	{
		iEdgeEnd += Abc_ObjFaninNum(pObj);
		if (iEdgeEnd - iEdgeBeg != 2)
		{
			Abc_Print(-1, "An internal AIG node is supposed to have 2 fanins.\n ");
			return 1;
		}
		iEdgeBeg = iEdgeEnd;
		iEdgeEnd += Abc_ObjFanoutNum(pObj);
		if (iEdgeEnd == iEdgeBeg)
		{
			Abc_Print(-1, "An internal AIG node is supposed to have at least 1 fanout.\n ");
			return 1;
		}
		pObj->iTemp = iNode++;
		m_vNodeIsPi.push_back(false);
		m_vpObjs.push_back(pObj);
		m_vpObjsNew.push_back(NULL);
		m_vNodeIndices.push_back(iEdgeEnd);
		iEdgeBeg = iEdgeEnd;

		pFanin0 = Abc_ObjFanin0(pObj);
		pFanin1 = Abc_ObjFanin1(pObj);
		
		if(fMetis)
		{
			weight = computeNodeWeightByLevel(pFanin0->Level, pFanin1->Level);
			//weight = computeNodeWeightByDegree(Abc_ObjFanoutNum(pObj));
			m_vNodeWeights.push_back(weight);
		}
	}
	Vec_PtrFree(vNodes);
	//COs
	Abc_NtkForEachCo(pNtk, pObj, i)
	{
		iEdgeEnd += Abc_ObjFaninNum(pObj);
		if (iEdgeEnd - iEdgeBeg != 1)
		{
			Abc_Print(-1, "A Co is supposed to have 1 fanin.\n ");
			return 1;
		}
		pObj->iTemp = iNode++;
		m_vNodeIsPi.push_back(false);
		m_vpObjs.push_back(pObj);
		m_vpObjsNew.push_back(NULL);
		m_vNodeIndices.push_back(iEdgeEnd);
		iEdgeBeg = iEdgeEnd;

		if(fMetis)
			m_vNodeWeights.push_back(1);
	}
	/*wxx
	Vec_PtrForEachEntry(Abc_Obj_t *, pMetis->AbcObj , pObj, i)
	{
		printf("i: %d\ttype: %d\n", i ,Abc_ObjType(pObj));
	}
	*/

	yassert(m_nNodes == m_vNodeIndices.size() - 1);
	yassert(m_nNodes == m_vNodeIsPi.size());
	yassert(m_nNodes == m_vpObjs.size());
	yassert(m_nNodes == m_vpObjsNew.size());
	yassert(m_nNodes == iNode);
//	yassert(m_nNodes == m_vNodeWeights.size());
	
	Abc_Print(-2, "\nThere are %d nodes in the initial graph to be partitioned.\n", m_nNodes);

	//add edges (m_vEdgeIsC is not prepared)
	for(i = 0; i < m_vpObjs.size(); i++)
	{
		pObj = m_vpObjs[i];
		if(isNodePi(i))
		{
			iEdge = m_vNodeIndices[i];
			Abc_ObjForEachFanout(pObj, pFanout, j)
			{
				m_vEdges.push_back(pFanout->iTemp);
				m_vEdgeIsC.push_back(false); //only first two edge for each non-Pi node has compl info
				iEdge++;

				if(fMetis)
					m_vEdgeWeights.push_back(1);
			}
			yassert(iEdge == m_vNodeIndices[i+1]);
		}
		else if(isNodeAND(i))
		{
			iEdge = m_vNodeIndices[i];
			pFanin0 = Abc_ObjFanin0(pObj);
			m_vEdges.push_back(pFanin0->iTemp);
			m_vEdgeIsC.push_back(pObj->fCompl0);

			if(fMetis)
			{
				weight = computeEdgeWeight(pFanin0->Level, nodeId2ReqTime[pFanin0->Id], pObj->Level, nodeId2ReqTime[pObj->Id]);
				m_vEdgeWeights.push_back(weight);
			}

			pFanin1 = Abc_ObjFanin1(pObj);
			m_vEdges.push_back(pFanin1->iTemp);
			m_vEdgeIsC.push_back(pObj->fCompl1);

			if(fMetis)
			{
				weight = computeEdgeWeight(pFanin1->Level, nodeId2ReqTime[pFanin1->Id], pObj->Level, nodeId2ReqTime[pObj->Id]);
				m_vEdgeWeights.push_back(weight);
			}

			iEdge += 2;
			Abc_ObjForEachFanout(pObj, pFanout, j)
			{
				m_vEdges.push_back(pFanout->iTemp);
				m_vEdgeIsC.push_back(false);
				if(fMetis)
				{
					weight = computeEdgeWeight(pObj->Level, nodeId2ReqTime[pObj->Id], pFanout->Level, nodeId2ReqTime[pFanout->Id]);
					m_vEdgeWeights.push_back(weight);
				}
				iEdge++;
			}
			yassert(iEdge == m_vNodeIndices[i+1]);
		}
		else if(isNodePo(i))
		{
			iEdge = m_vNodeIndices[i];
			pFanin0 = Abc_ObjFanin0(pObj);
			m_vEdges.push_back(pFanin0->iTemp);
			m_vEdgeIsC.push_back(pObj->fCompl0);
			iEdge++;
			//weight = computeEdgeWeight(pFanin0->Level, nodeId2ReqTime[pFanin0->Id], pObj->Level, nodeId2ReqTime[pObj->Id]);
			if(fMetis)
				m_vEdgeWeights.push_back(WALL_EDGE_WEIGHT);
			yassert(iEdge == m_vNodeIndices[i+1]);
		}
		else
			yassert(0);
	}
	yassert(m_vEdges.size() == m_vEdgeIsC.size());
	//yassert(m_vEdges.size() == m_vEdgeWeights.size());

	ylog("Total number of edges in MetisGraph: %ld\n", m_vEdges.size()/2);
	return 0;
}

/* comment metis by zli
int MetisGraph::partGraphByMetis(int32_t nParts)
{
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options); 
	options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
	options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
	//options[METIS_OPTION_NO2HOP] = 1;
	options[METIS_OPTION_UFACTOR] = METIS_UFACTOR;
	m_nParts = nParts;

	struct timeval t1,t2;
	double time;
	gettimeofday(&t1, NULL);

	if(m_vEdgeWeights.empty() || m_vNodeWeights.empty() )
		METIS_PartGraphKway(&m_nNodes, &m_nCon, &m_vNodeIndices[0], &m_vEdges[0], NULL, NULL, NULL, &m_nParts, NULL, NULL, NULL, &m_nCutEdges, &m_vPartition[0]);
	else
		METIS_PartGraphKway(&m_nNodes, &m_nCon, &m_vNodeIndices[0], &m_vEdges[0], &m_vNodeWeights[0], NULL, &m_vEdgeWeights[0], &m_nParts, NULL, NULL, options, &m_nCutEdges, &m_vPartition[0]);

	gettimeofday(&t2, NULL);
   	time = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
   	printf("METIS spent time: %f\n", time);
	Abc_Print(-2, "nedgecut = %d\n", m_nCutEdges);
	return m_nCutEdges;
}
*/

/**Function*************************************************************

  Synopsis    [Start a subnetwork using exsiting network.]

  Description [Duplicate CI/COs from original network, leaving out internal nodes and newly generated PI/POs.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
Abc_Ntk_t* MetisGraph::initOneSubNtk(int index)
{
	Abc_Ntk_t *pNtkNew;
	Abc_Obj_t *pObj;
	int i, ipart;
	// int fCopyNames;

	if (m_pNtkOrigin == NULL)
		return NULL;
	yassert(index >= 0 && index < m_nParts);
	// decide whether to copy the names
	// fCopyNames = (pNtk->ntkType != ABC_NTK_NETLIST);
	
	// start the network
	pNtkNew = Abc_NtkAlloc(m_pNtkOrigin->ntkType, m_pNtkOrigin->ntkFunc, 1);
	pNtkNew->nConstrs = m_pNtkOrigin->nConstrs;
	pNtkNew->nBarBufs = m_pNtkOrigin->nBarBufs;

	//zkx
	pNtkNew->ntkType = m_pNtkOrigin->ntkType;

	//
	// duplicate the name and the spec
	pNtkNew->pName = Extra_UtilStrsav(m_pNtkOrigin->pName);
	pNtkNew->pSpec = Extra_UtilStrsav(m_pNtkOrigin->pSpec);
	// clone CI/CO Nodes
	for(i = 0; i < m_nNodes; i++)
	{
		ipart = m_vPartition[i];
		if (ipart != index)
			continue;
		pObj = m_vpObjs[i];
		if (Abc_ObjType(pObj) == ABC_OBJ_NODE)
			continue;
		//yassert( pObj->iTemp == i );
		if (Abc_ObjType(pObj) == ABC_OBJ_CONST1)
			pObj->pCopy = Abc_AigConst1(pNtkNew);
		else if ( Abc_ObjIsCi(pObj) )
		{
			pObj->pCopy = Abc_NtkCreatePi(pNtkNew);
			Abc_ObjAssignName( pObj->pCopy, Abc_ObjName(pObj), NULL);
			yassert(pObj->Level == 0);
			pObj->pCopy->Level = pObj->Level;
		}
		else if ( Abc_ObjIsCo(pObj))
		{
			pObj->pCopy = Abc_NtkCreatePo(pNtkNew);
			Abc_ObjAssignName( pObj->pCopy, Abc_ObjName(pObj), NULL);				
		}
		else 
			yassert(0);
		pObj->pCopy->fCompl0 = pObj->fCompl0;
		pObj->pCopy->fCompl1 = pObj->fCompl1;
		//store the newly created AbcObj for corresponding node in Metis graph
		m_vpObjsNew[i] = pObj->pCopy;
	}
	Abc_ManTimeDup( m_pNtkOrigin, pNtkNew );
    pNtkNew->AndGateDelay = m_pNtkOrigin->AndGateDelay;

    if ( pNtkNew->pManTime && Abc_FrameReadLibGen() && pNtkNew->AndGateDelay == 0.0 )
	{
		pNtkNew->AndGateDelay = Mio_LibraryReadDelayAigNode((Mio_Library_t *)Abc_FrameReadLibGen());  
	}	  
	return pNtkNew;
}

/**Function*************************************************************

  Synopsis    [Generate subnetworks using hmetis partition results.]

  Description [Create internal AND nodes and new PI/PO pairs for edges that straddle different partitions in DFS order(original pMetis->AbcObj order). Finish the connections as well.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
int MetisGraph::createSubNtksFromPartition(vector<Abc_Ntk_t*>& vSubNtks)
{
	Abc_Obj_t *pObj, *pFanin0, *pFanin1, *pFanout;
	Abc_Obj_t *pObjNew, *pPiNew, *pPoNew;
	Abc_Ntk_t *pNtk;
	int iFanin0, iFanin1, ipart, ipart0, ipart1;
	int i, j;
	cout << "vSubNtks.size() is " << vSubNtks.size() << " , m_nParts is " << m_nParts << endl;


	yassert(vSubNtks.size() == m_nParts);
	int nCutEdgeFromPi = 0;
	int nCutEdge= 0;

	for(i = 0; i < vSubNtks.size(); i++)
	{
		vSubNtks[i] = initOneSubNtk(i);
	}

	//DFS order is expected.
	for(i = 0; i < m_vpObjs.size(); i++)
	{
		pObj = m_vpObjs[i];
		if (Abc_ObjIsCi(pObj) || pObj->Type == ABC_OBJ_CONST1)
			continue;
		ipart = m_vPartition[i];
		pNtk = vSubNtks[ipart];//the pSubNtk that the current pObj belongs to

		//check fanin0 for CO and And nodes
		iFanin0 = getNodeFanin0Id(i);
		pFanin0 = m_vpObjsNew[iFanin0]; //pObj in pSubNtk
		if (pFanin0 == NULL)
		{
			Abc_Print(-1, "The AIG Nodes in metis array are expected to be in DFS order.\n ");
			return 1;
		}
		ipart0 = m_vPartition[iFanin0];
		if (ipart0 != ipart) //cut edge found!
		{
			nCutEdge++;
			if(Abc_ObjIsPi(pFanin0)) //when the Ci has a fanout edge that is cut, don't generate cut-caused PO
			{
				nCutEdgeFromPi++;
				pPiNew = Abc_NtkCreatePi(pNtk);
				Abc_ObjAssignName(pPiNew, Abc_ObjName(pPiNew), NULL);
				pPiNew->fMarkA = 1;
				pPiNew->pData = pFanin0; //map the cut-caused PI to the corresponding PO
				pFanin0 = pPiNew;		//point the pFanin0 to the corresponding fanin pObj in the same pSubNtk
			}
			else
			{
				pPoNew = NULL;
				Abc_ObjForEachFanout(pFanin0, pFanout, j)
				{
					//If cut-caused PO already exists, use it directly, avoiding duplication.
					if (Abc_ObjIsPo(pFanout) && pFanout->fMarkA == 1)
					{
						pPoNew = pFanout;
						break;
					}
				}
				if (pPoNew == NULL) //create a new cut-caused PO
				{
					pPoNew = Abc_NtkCreatePo(vSubNtks[ipart0]);
					Abc_ObjAssignName(pPoNew, Abc_ObjName(pPoNew), NULL); //assign a temporary name
					pPoNew->fMarkA = 1;									  //mark as cut-caused PO
					Abc_ObjAddFanin(pPoNew, pFanin0);
				}
				pPiNew = Abc_NtkCreatePi(pNtk);
				Abc_ObjAssignName(pPiNew, Abc_ObjName(pPiNew), NULL);
				pPiNew->fMarkA = 1;
				pPiNew->pData = pPoNew; //map the cut-caused PI to the corresponding PO
				pFanin0 = pPiNew;		//point the pFanin0 to the corresponding fanin pObj in the same pSubNtk
			}

		}
		if (Abc_ObjIsCo(pObj))
		{
			pObjNew = m_vpObjsNew[i]; //pObjNew of CI/CO in pSubNtk has been added in Abc_SubNtkStartFrom()
			yassert(pObjNew != NULL);
			Abc_ObjAddFanin(pObjNew, pFanin0); //Connect
		}
		else if (Abc_ObjIsNode(pObj))
		{
			//fanin1
			iFanin1 = getNodeFanin1Id(i);
			pFanin1 = m_vpObjsNew[iFanin1]; //pObj in pSubNtk. In DFS, fanin should have been processed
			yassert(pFanin1);
			if (pFanin1 == NULL)
			{
				Abc_Print(-1, "The AIG Nodes in metis array are expected to be in DFS order.\n ");
				return 1;
			}
			ipart1 = m_vPartition[iFanin1];
			if (ipart1 != ipart)
			{
				nCutEdge++;
				if(Abc_ObjIsPi(pFanin1)) //when the Ci has a fanout edge that is cut, don't generate cut-caused PO
				{
					nCutEdgeFromPi++;
					pPiNew = Abc_NtkCreatePi(pNtk);
					Abc_ObjAssignName(pPiNew, Abc_ObjName(pPiNew), NULL);
					pPiNew->fMarkA = 1;
					pPiNew->pData = pFanin1; //map the cut-caused PI to the corresponding PO
					pFanin1 = pPiNew;		//point the pFanin0 to the corresponding fanin pObj in the same pSubNtk
				}
				else
				{
					pPoNew = NULL;
					Abc_ObjForEachFanout(pFanin1, pFanout, j)
					{
						if (Abc_ObjIsPo(pFanout) && pFanout->fMarkA == 1)
						{
							pPoNew = pFanout;
							break;
						}
					}
					if (pPoNew == NULL)
					{
						pPoNew = Abc_NtkCreatePo(vSubNtks[ipart1]);
						Abc_ObjAssignName(pPoNew, Abc_ObjName(pPoNew), NULL);
						pPoNew->fMarkA = 1;
						Abc_ObjAddFanin(pPoNew, pFanin1);
					}
					pPiNew = Abc_NtkCreatePi(pNtk);
					Abc_ObjAssignName(pPiNew, Abc_ObjName(pPiNew), NULL);
					pPiNew->fMarkA = 1;
					pPiNew->pData = pPoNew;
					pFanin1 = pPiNew;
				}
			}

			pObjNew = Abc_AigAnd(static_cast<Abc_Aig_t*>(pNtk->pManFunc), Abc_ObjNotCond(pFanin0, isNodeFanin0C(i)), Abc_ObjNotCond(pFanin1, isNodeFanin1C(i)));

			m_vpObjsNew[i] = pObjNew;
			pObj->pCopy = pObjNew;
		}
		else
			yassert(0);
	}

	
	ylog("Total number of CIs with cut fanout edges: %d\n", nCutEdgeFromPi);
	ylog("Total number of cut edges: %d\n", nCutEdge);

	return 0;
}


void MetisGraph::tmp()
{
	Yaig y;
	struct timeval t1,t2;
	double time;
	gettimeofday(&t1, NULL);
	y.addSubGraph(*this);
	gettimeofday(&t2, NULL);
   	time = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
   	printf("yaig::addSubGraph spent time: %f\n", time);
}

void MetisGraph::setPoPart()
{
	for(int i = 0; i < m_nNodes; i++)
	{
		if(m_vPartition[i] == -1)
		{
			if(!isNodePo(i))
			{
				if(isNodePi(i))
					ylog("Pi has no partition!!\n");
				else
					ylog("AND node[%d] has no partition!!\n", i);
				yassert(0);
			}
			int faninId = getNodeFanin0Id(i);
			yassert(m_vPartition[faninId] != -1);
			m_vPartition[i] = m_vPartition[faninId];
		}
	}
	return;
}

void MetisAig::bindGraph(MetisGraph* pmg)
{
	if(m_pMG)
	{
		clear();
		vector<Cone>().swap(m_vCones);
		vector<Cluster>().swap(m_vClusters);
		vector<int>().swap(m_vConeId2ClusterId);
		m_pMG = NULL;
	}
	addSubGraph(*pmg);
	m_pMG = pmg;
	return ;
}

void MetisAig::printCones()
{
	ylog("There are %ld cones\n", m_vCones.size());
 #if 0
	for(auto& c : m_vCones)
	{
		if(c.vBoundaryEdges.size())
			ylog("m_vCones[%d] has %ld boundary edges, level: %d\n",c.iId, c.vBoundaryEdges.size(), c.iMaxLevel);
		for(auto e : c.vBoundaryEdges)
			ylog("edge: %d -> %d\n", e.iFaninId, e.iFanoutId);
		printOneCone(c.iId);
	}
 #endif
}
 
void MetisAig::printOneCone(int32_t coneId)
{
	ylog("Now print m_vCones[%d]:\n", coneId);
	auto& cone = m_vCones[coneId];
	ylog("iMaxLevel = %d, iPoId = %d, iId = %d\n", cone.iMaxLevel, cone.iPoId, cone.iId);
	ylog("NumBoundaryEdges = %ld, NumNeighborCone = %ld\n", cone.vBoundaryEdges.size(), cone.sAdjacentConeId.size());
}

void MetisAig::printClusters()
{
	ylog("There are %ld clusters\n", m_vClusters.size());
	for(int i = 0; i <  m_vClusters.size(); i++)
	//std::min(10, m_vClusters.size())
		ylog("cluster[%d]: workload = %d\tiMaxLevel = %d\tnNodes = %d\tnCones = %ld\n", i, m_vClusters[i].iWorkload, m_vClusters[i].iMaxLevel, m_vClusters[i].nNodes, m_vClusters[i].vConeIds.size()) ;
 #if 0 
	for(auto& c : m_vClusters)
	{
		if(c.vConeIds.size())
			ylog("this cluster has %ld cones\n", c.vConeIds.size());
		for(auto coneId : c.vConeIds)
			ylog("coneId: %d\n", coneId);
	}
 #endif
}

//将AIG电路的节点划分成多个cluster
void MetisAig::parseAig()
{
	//yassert(m_pMG);
	//prepare m_vCones
	//计算所有节点的level
	computeAllLevel();

	//对于每一个PO节点，构建一个Cone
	for(auto poId : m_vPos)
	{
		Cone cone(m_vNodes[poId].iLevel, poId);
		m_vCones.push_back(move(cone));
	}

	//排序，按照level大小
	sort(m_vCones.begin(), m_vCones.end(), [](const Cone& c1, const Cone& c2){ return c1.iMaxLevel > c2.iMaxLevel;});
	for(int i = 0; i < m_vCones.size(); i++)
	{
		auto& cone = m_vCones[i];
		cone.iId = i;
		//找到它的边界边（即与其他Cone对象相邻的边界边），构建一个BoundaryEdge对象
		findBoundaryEdges(cone);
	}
	//printCones();

	for(auto& cone : m_vCones)
		m_tmp0 += cone.vBoundaryEdges.size();
	ylog("There are totally %d boundary edges\n", m_tmp0);
	m_tmp0 = 0;

	vector<int> coneId2ClusterId(m_vCones.size(), -1);

	//prepare m_vClusters
	//根据每个Cone对象的相邻关系，将它们划分到多个Cluster对象中
	for(auto& cone : m_vCones)
	{
		if(m_vClusters.empty())
		{
			Cluster cluster(0, cone.iMaxLevel);
			cluster.vConeIds.push_back(cone.iId);
			coneId2ClusterId[cone.iId] = cluster.iId;
			m_vClusters.push_back(move(cluster));
			continue;
		}
		for(auto neighborConeId : cone.sAdjacentConeId)
		{
			int32_t clusterIdNeighborCone = coneId2ClusterId[neighborConeId];
			yassert(clusterIdNeighborCone != -1);
			int32_t clusterIdCurrentCone = coneId2ClusterId[cone.iId];
			if(clusterIdCurrentCone == -1)
			{
				coneId2ClusterId[cone.iId] = clusterIdNeighborCone; //longer cone should belong to existed cluster
				Cluster& cluster = m_vClusters[clusterIdNeighborCone];
				cluster.vConeIds.push_back(cone.iId);
			}
			else if(clusterIdCurrentCone != clusterIdNeighborCone)
			{
				Cluster& cluster = m_vClusters[clusterIdNeighborCone];
				Cluster& clusterCur = m_vClusters[clusterIdCurrentCone];
				//clusterCur merges to cluster
				clusterCur.iId = -1;
				if(cluster.iMaxLevel < clusterCur.iMaxLevel)
					cluster.iMaxLevel = clusterCur.iMaxLevel;
				for(auto coneId : clusterCur.vConeIds)
				{
					yassert(coneId2ClusterId[coneId] == clusterIdCurrentCone);
					coneId2ClusterId[coneId] = clusterIdNeighborCone;
					cluster.vConeIds.push_back(coneId);
				}
			}
		}
		if(coneId2ClusterId[cone.iId] == -1)
		{
			Cluster cluster(m_vClusters.size(), cone.iMaxLevel);
			cluster.vConeIds.push_back(cone.iId);
			coneId2ClusterId[cone.iId] = cluster.iId;
			m_vClusters.push_back(move(cluster));
		}
	}
	//printClusters();

	//compute clusters' workload
	m_iTotalWorkLoad = 0;
	m_iMaxClusterWorkLoad = 0;
	for(auto& cluster : m_vClusters)
	{
		if(cluster.iId == -1) //merged cluster
		{
			cluster.iWorkload = -1;
			cluster.nNodes = -1;
			continue;
		}
		computeWorkLoad(cluster);
		m_iTotalWorkLoad += cluster.iWorkload;
	}
	//sort clusters by workload
	//after sorting, the iId field of each cluster in the m_vClusters should be re-assigned. 
	sort(m_vClusters.begin(), m_vClusters.end(), [](const Cluster& lhs, const Cluster& rhs){return lhs.iWorkload > rhs.iWorkload;});
	m_iMaxClusterWorkLoad = m_vClusters[0].iWorkload;
	//m_vConeId2ClusterId.resize(m_vCones.size(), -1);
	for(int i = 0; i < m_vClusters.size(); i++)
	{
		if(m_vClusters[i].iId == -1)
		{
			m_vClusters.resize(i);
			break;
		}
		auto& cluster = m_vClusters[i];
		cluster.iId = i;
#if 0
		for(auto coneId : cluster.vConeIds)
		{
			yassert(m_vConeId2ClusterId[coneId] == -1);
			m_vConeId2ClusterId[coneId] = cluster.iId;
		}
#endif
	}
	
	// merge small cluster
	uint32_t ClusterUpB = m_pMG->get_sCluster();
	if (ClusterUpB != 0)
		mergeSmallClusters(ClusterUpB);
	
	printClusters();
#if 0 
	//checkClusters();
	//for(int i = 0; i < m_vClusters.size() ; i++)
	for(int i = 0; i < 20 ; i++)
		ylog("cluster[%d]: workload = %d\tiMaxLevel = %d\tnNodes = %d\tnCones = %ld\n", i, m_vClusters[i].iWorkload, m_vClusters[i].iMaxLevel, m_vClusters[i].nNodes, m_vClusters[i].vConeIds.size());
#endif
	ylog("total work load: %d\n", m_iTotalWorkLoad);
	ylog("max cluster work load: %d\n", m_iMaxClusterWorkLoad);
}

void MetisAig::mergeSmallClusters(uint32_t size) {
	sort(m_vClusters.begin(), m_vClusters.end(), [](const Cluster& lhs, const Cluster& rhs){return lhs.nNodes < rhs.nNodes;});
	vector<Cluster> merged_vSmallCluters;
	for (auto& cluster : m_vClusters) {
		if (merged_vSmallCluters.empty())
			merged_vSmallCluters.push_back(cluster);
		else {
			if (merged_vSmallCluters.back().nNodes + cluster.nNodes > size)
				merged_vSmallCluters.push_back(cluster);
			else {
				merged_vSmallCluters.back().iMaxLevel = max(merged_vSmallCluters.back().iMaxLevel, cluster.iMaxLevel);
				merged_vSmallCluters.back().iWorkload += cluster.iWorkload;
				merged_vSmallCluters.back().nNodes += cluster.nNodes;
				merged_vSmallCluters.back().vConeIds.insert(merged_vSmallCluters.back().vConeIds.end(), cluster.vConeIds.begin(), cluster.vConeIds.end());
			}
		}
	}
	m_vClusters = merged_vSmallCluters;
	sort(m_vClusters.begin(), m_vClusters.end(), [](const Cluster& lhs, const Cluster& rhs){return lhs.iWorkload > rhs.iWorkload;});
	m_iMaxClusterWorkLoad = m_vClusters[0].iWorkload;
	for(int i = 0; i < m_vClusters.size(); i++)
		m_vClusters[i].iId = i;
}

void MetisAig::checkClusters()
{
	for(auto& cluster : m_vClusters)
	{
		setNextIter();
		for(auto coneId : cluster.vConeIds)
		{
			auto& cone = m_vCones[coneId];
			checkNodeCluster(cone.iPoId, cluster.iId);
		}
	}
}

void MetisAig::checkNodeCluster(int32_t nodeId, int32_t clusterId)
{
	Node& node = m_vNodes[nodeId];
	if(node.isPo())
	{
		checkNodeCluster(node.getFanin0Id(), clusterId);
		return;
	}
	if(node.isPi())
		return;
	if(node.iIter == m_iGlobalIter)
		return;
	node.iIter = m_iGlobalIter;
	yassert(m_vConeId2ClusterId[node.iConeId] == clusterId);
	checkNodeCluster(node.getFanin0Id(), clusterId);
	checkNodeCluster(node.getFanin1Id(), clusterId);
	return;
}


void MetisAig::findBoundaryEdges(Cone& cone)
{
	setNextIter();
	visitAllFaninFromNode(cone.iPoId, cone); //mark all the AND nodes in this cone and prepare nVisits of the nodes
	setNextIter();
	findBoundaryEdges_rec(cone, cone.iPoId, 0);
	return;
}

//遍历以指定节点为起点的所有 fanin 节点，同时更新每个节点的访问状态和节点所属的cone区域
void MetisAig::visitAllFaninFromNode(int32_t nodeId, Cone& cone)
{
	Node& node = m_vNodes[nodeId];
	if(node.isPo())
	{
		visitAllFaninFromNode(node.getFanin0Id(), cone);
		return;
	}
	if(node.iIter == m_iGlobalIter)
	{
		yassert(node.nVisits);
		node.nVisits--;
		return ;
	}
	/*
	if(node.isPi())
		return ;
	*/
	node.iIter = m_iGlobalIter;
	node.nVisits = node.nFanouts - 1;
	if(node.iConeId == -1)
		node.iConeId = cone.iId;
	else
		if(cone.iMaxLevel > m_vCones[node.iConeId].iMaxLevel)
			yassert(0);
	if(node.isPi())
		return;
			//node.iConeId = cone.iId;
	visitAllFaninFromNode(node.getFanin0Id(), cone);
	visitAllFaninFromNode(node.getFanin1Id(), cone);
	return;
}

//在一个有向无环图中找到一个锥形子图的边界边
void MetisAig::findBoundaryEdges_rec(Cone& cone, int32_t nodeId, int32_t fCovered)
{
	Node& node = m_vNodes[nodeId];
	int fc = fCovered;

	//检查是否为PO
	if(node.isPo())
	{
		//递归访问它的输入节点，并检查是否为边界边
		Node& fanin = m_vNodes[node.getFanin0Id()];
		if(fanin.nVisits)// && !fanin.isPi())
		{
			Edge edge(fanin.id, nodeId);
			cone.vBoundaryEdges.push_back(move(edge));
			if(fanin.iConeId != cone.iId)
				cone.sAdjacentConeId.insert(fanin.iConeId);
			fc = 1;
		}
		findBoundaryEdges_rec(cone, node.getFanin0Id(), fc);
		return;
	}
	if(node.isPi())
		return;
	if(node.iIter == m_iGlobalIter)
		return;
	node.iIter = m_iGlobalIter;

	//如果不是输出节点，那么它会继续访问两个输入节点，并对它们进行递归调用。
	Node& fanin0 = m_vNodes[node.getFanin0Id()];
	Node& fanin1 = m_vNodes[node.getFanin1Id()];

	//如果该节点的一个输入节点已被访问过或它的多个输出节点不属于给定锥体(cone)，那么该节点就是边界节点。
	if(fanin0.nVisits || (fanin0.nFanouts > 1 && fanin0.iConeId != cone.iId))// && !fanin0.isPi()) //has fanout node that belongs other cone
	{
		if(!fCovered)
		{
			//每个边界边被保存到锥体的边界边列表中
			Edge edge(fanin0.id, nodeId);
			cone.vBoundaryEdges.push_back(move(edge));
		}
		if(fanin0.iConeId != cone.iId)
			//与其相邻的其他锥体的ID被保存到锥体的相邻锥体ID集合中
			cone.sAdjacentConeId.insert(fanin0.iConeId);
		fc = 1;
	}
	findBoundaryEdges_rec(cone, fanin0.id, fc);

	fc = fCovered;
	if(fanin1.nVisits || (fanin1.nFanouts > 1 && fanin1.iConeId != cone.iId))// && !fanin1.isPi()) //has fanout node that belongs other cone
	{
		if(!fCovered)
		{
			Edge edge(fanin1.id, nodeId);
			cone.vBoundaryEdges.push_back(move(edge));
		}
		if(fanin1.iConeId != cone.iId)
			cone.sAdjacentConeId.insert(fanin1.iConeId);
		fc = 1;
	}
	findBoundaryEdges_rec(cone, fanin1.id, fc);

	return;
}

void MetisAig::computeWorkLoad(Cluster& cluster)
{
	cluster.iWorkload = 0;
	cluster.nNodes = 0;
	setNextIter();
	for(auto coneId : cluster.vConeIds)
		computeWorkLoad_rec(m_vCones[coneId].iPoId, cluster.iWorkload, cluster.nNodes);
}

int32_t MetisAig::computeWorkLoad_rec(int32_t nodeId, int32_t& workload, int32_t& nNodes)
{
	Node& node = m_vNodes[nodeId];
	if(node.isPo())
	{
		computeWorkLoad_rec(node.getFanin0Id(), workload, nNodes);
		return workload;
	}
	if(node.isPi())
	{
		workload += PI_WORK_LOAD;
		return 1;
	}
	if(node.iIter == m_iGlobalIter)
		return node.iNCuts;
	node.iIter = m_iGlobalIter;
	nNodes++;

	int32_t nCuts, nCutsLeft, nCutsRight;
	nCutsLeft = computeWorkLoad_rec(node.getFanin0Id(), workload, nNodes);
	nCutsRight = computeWorkLoad_rec(node.getFanin1Id(), workload, nNodes);
	nCuts = nCutsLeft * nCutsRight;
	workload += nCuts;
	if(nCuts > 8)
		nCuts = 8;
	node.iNCuts = nCuts + 1;
	return node.iNCuts;
}

int32_t MetisAig::decideNumParts()
{
	int nParts = m_iTotalWorkLoad / m_iMaxClusterWorkLoad;
	ylog("Max nParts = %d\n", nParts);
	if(nParts < MIN_N_PART)
	{
		//nParts = tryPart();
		//nParts = tryPart2();
		ylog("This graph cannot be partitioned naturally. Metis will be used\n");
		return -1;
	}
	else 
	{
		while(nParts > MAX_N_PART)
			nParts = nParts / 2;
	}
	ylog("At last, nParts = %d\n", nParts);
	return nParts;
}

int32_t MetisAig::tryPart2()
{
	int clusterId = 0;
	int newClusterId = m_vClusters.size();
	sort(m_vClusters[clusterId].vConeIds.begin(), m_vClusters[clusterId].vConeIds.end(), [this](int32_t lhs, int32_t rhs){
		return m_vCones[lhs].iId < m_vCones[rhs].iId;
	});
	for(int i = 0; i < 100; i++)
		printOneCone(m_vClusters[clusterId].vConeIds[i]);

	vector<int>(m_vCones.size(), -1).swap(m_vConeId2ClusterId);
	for(auto coneId : m_vClusters[clusterId].vConeIds)
		m_vConeId2ClusterId[coneId] = m_vClusters[clusterId].iId;
	vector<int> coneIdForCut;

	auto& coneCritical = m_vCones[m_vClusters[clusterId].vConeIds[0]];
	coneIdForCut.push_back(coneCritical.iId);
	m_vConeId2ClusterId[coneCritical.iId] = newClusterId;
	//Not finished. This method is abandoned.

	return 2;

}

int32_t MetisAig::tryPart()
{
	int workLoadLimit = m_iTotalWorkLoad / MIN_N_PART;
	int iCluster = 0;
	int nConesCut;
	int nCluster = m_vClusters.size();
	//try PI-Cut
	while(m_iMaxClusterWorkLoad > workLoadLimit) 
	{
		//Cluster& cluster = m_vClusters[iCluster];
		nConesCut = partBiggestClusterByPICut(iCluster);
		if(nConesCut)
		{
			//if the next cluster becomes the biggest cluster:
			if(m_iMaxClusterWorkLoad == m_vClusters[++iCluster].iWorkload)
				workLoadLimit = m_iTotalWorkLoad / MIN_N_PART;
			else //if after Pi-Cut, the current cluster is still the biggest cluster:
				break;
		}
		else //The biggest cluster cannot be cut by PI-Cut
		{
			ylog("No Pi-Cut available for cluster[%d]\n", iCluster);
			break;
		}
	}
	if(nCluster < m_vClusters.size()) //if m_Clusters is changed, update the order and iId field.
	{
		sort(m_vClusters.begin(), m_vClusters.end(), [](const Cluster& lhs, const Cluster& rhs){return lhs.iWorkload > rhs.iWorkload;});
		for(int i = 0; i < m_vClusters.size(); i++)
			m_vClusters[i].iId = i;
		int totalwl = 0;
		for(auto& cluster : m_vClusters)
			totalwl += cluster.iWorkload;
		yassert(totalwl == m_iTotalWorkLoad);
		ylog("totalwl = %d\n", totalwl);
	}
	yassert(m_iMaxClusterWorkLoad == m_vClusters[0].iWorkload);
	printClusters();

	//return 2;

	iCluster = 0;
	//try Edge-Cut
	while(m_iMaxClusterWorkLoad > workLoadLimit)
	{
		nConesCut = partBiggestCluster(iCluster, workLoadLimit);
		printClusters();
		sort(m_vClusters.begin(), m_vClusters.end(), [](const Cluster& lhs, const Cluster& rhs){return lhs.iWorkload > rhs.iWorkload;});
		m_vConeId2ClusterId.resize(m_vCones.size(), -1);
		for(int i = 0; i < m_vClusters.size(); i++)
		{
			Cluster& cluster = m_vClusters[i];
			m_vClusters[i].iId = i;
			for(auto coneId : cluster.vConeIds)
			{
				yassert(m_vConeId2ClusterId[coneId] == -1);
				m_vConeId2ClusterId[coneId] = cluster.iId;
			}
		}
		checkClusters();

		return 2;
		yassert(0);
		break;
	}
	return 2;

}



int32_t MetisAig::partBiggestCluster(int32_t clusterId, int32_t workLoadLimit)
{
	if(m_vClusters[clusterId].iWorkload <= workLoadLimit)
		return 0;
	yassert(m_vClusters[clusterId].iWorkload == m_iMaxClusterWorkLoad);
	sort(m_vClusters[clusterId].vConeIds.begin(), m_vClusters[clusterId].vConeIds.end(), [this](int32_t lhs, int32_t rhs){
		return m_vCones[lhs].sAdjacentConeId.size() > m_vCones[rhs].sAdjacentConeId.size();
	});
	int nConesCut = 0;
	int workloadNew = 0;
	int workloadOld = m_vClusters[clusterId].iWorkload;
	int workloadNewMax = 0;
	for(int i = 0; i < m_vClusters[clusterId].vConeIds.size(); i++)
	{
		auto& cone = m_vCones[m_vClusters[clusterId].vConeIds[i]];
		bool fCut = 1;
		for(auto edge : cone.vBoundaryEdges)
		{
			auto& fanin = m_vNodes[edge.iFaninId];
			auto& fanout = m_vNodes[edge.iFanoutId];
#if 0
			if(fanout.iLevel - fanin.iLevel > 129)
			{
				if(fanout.iLevel - fanin.iLevel == 130)
				{
					ylog("===========\nymc: This may be a critical cut\n");
					ylog("cone[%d]: edge: node[%d]->node[%d]\n", cone.iId, fanin.id, fanout.id);
					m_tmp0++;
				}
				
				fCut = 0;
				break;
			}
#endif
			if(fanin.iConeId == cone.iId)
			{
				fCut = 0;
				break;
			}
			if(fanout.isPo())
			{
				fCut = 0;
				break;
			}
		}
		if(fCut) //Cut this cone from cluster
		{
			for(auto edge : cone.vBoundaryEdges)
				cutOneBoundaryEdge(edge);
			//vector<Edge>().swap(cone.vBoundaryEdges);
			//set<int32_t>().swap(cone.sAdjacentConeId);
			cone.iMaxLevel = computeNodeLevel(m_vNodes[cone.iPoId].getFanin0Id());

			Cluster clusterNew(m_vClusters.size(), cone.iMaxLevel);
			clusterNew.vConeIds.push_back(cone.iId);
			computeWorkLoad(clusterNew);
			if(workloadNewMax < clusterNew.iWorkload)
				workloadNewMax = clusterNew.iWorkload;
			workloadNew += clusterNew.iWorkload;
			m_vClusters.push_back(move(clusterNew));
			nConesCut++;

			m_vClusters[clusterId].vConeIds[i] = -1;	
		}
	}

	if(nConesCut) //EdgeCut happened!
	{
		ylog("In partCluster(): %d EdgeCut happened for cluster[%d]\n", nConesCut, clusterId);
		int maxLev = 0;
		sort(m_vClusters[clusterId].vConeIds.begin(), m_vClusters[clusterId].vConeIds.end(), [](int32_t lhs, int32_t rhs){ return lhs > rhs; });
		for(int i = 0; i < m_vClusters[clusterId].vConeIds.size(); i++)
		{
			if(m_vClusters[clusterId].vConeIds[i] == -1) //delete cut cones
			{
				m_vClusters[clusterId].vConeIds.resize(i);
				break;
			}
			auto& cone = m_vCones[m_vClusters[clusterId].vConeIds[i]];
			if(maxLev < cone.iMaxLevel) //update maxlevel
				maxLev = cone.iMaxLevel;
		}
		if(m_vClusters[clusterId].iMaxLevel > maxLev) //critical cone is parted away
			m_vClusters[clusterId].iMaxLevel = maxLev;
		computeWorkLoad(m_vClusters[clusterId]);
		m_iTotalWorkLoad = m_iTotalWorkLoad - workloadOld + workloadNew + m_vClusters[clusterId].iWorkload;
		m_iMaxClusterWorkLoad = m_vClusters[clusterId].iWorkload; //assume the current cluster is still the biggest one.
		if(m_iMaxClusterWorkLoad <  workloadNewMax) //cut-out cone may be the biggest cluster?
			m_iMaxClusterWorkLoad = workloadNewMax;
		if(m_iMaxClusterWorkLoad < m_vClusters[clusterId + 1].iWorkload) //After PI-Cut, next cluster may be bigger?
			m_iMaxClusterWorkLoad = m_vClusters[clusterId + 1].iWorkload; 
	}

	if(!nConesCut)
		ylog("In partCluster(): No EdgeCut happened for cluster[%d]\n", clusterId);

	ylog("After partCluster(), m_iTotalwl = %d, m_iMaxwl = %d\n", m_iTotalWorkLoad, m_iMaxClusterWorkLoad);

	if(m_iMaxClusterWorkLoad > workLoadLimit)
		ylog("After partCluster(), the workLoadLimit can not be satisfied.\n");

	return nConesCut;
}

void MetisAig::cutOneBoundaryEdge(Edge& edge)
{
	//Warning: this function may destroy the DFS order of m_vNodes, and invalidate the structure hash table.
	//The connectivity of nodes is modified directly
	Node& fanin = m_vNodes[edge.iFaninId];
	Node& fanout = m_vNodes[edge.iFanoutId];
	int32_t newPiId = addPi();
	Literal fanin0, fanin1;
	if(fanout.getFanin0Id() == fanin.id)
		fanin0 = Literal(newPiId, fanout.isC0());
	else
		fanin0 = fanout.attr.getFanin0();

	if(fanout.getFanin1Id() == fanin.id)
		fanin1 = Literal(newPiId, fanout.isC1());
	else
		fanin1 = fanout.attr.getFanin1();
	
	NodeAttribute nodeAttr(fanin0, fanin1);
	fanout.attr = nodeAttr;
}

int32_t MetisAig::partBiggestClusterByPICut(int32_t clusterId)
{
	yassert(m_vClusters[clusterId].iWorkload == m_iMaxClusterWorkLoad);
	yassert(m_vClusters[clusterId].vConeIds.size() > 1);
	int32_t res = 0; //the number of the cut cones
	int32_t workloadNewMax = 0;
	bool fPiCut;
	for(auto i = 0; i < m_vClusters[clusterId].vConeIds.size(); i++)
	{
		auto& cone = m_vCones[m_vClusters[clusterId].vConeIds[i]];
		yassert(!cone.vBoundaryEdges.empty());
		fPiCut = true;
		for(auto edge : cone.vBoundaryEdges) //make sure all the boundary edges of the cone are fanouts of PI
		{
			auto fanin = m_vNodes[edge.iFaninId];
			if(!fanin.isPi())
			{
				fPiCut = false;
				break;
			}
		}
		if(fPiCut)
		{
			//generate new cluster from the cut cone
			Cluster clusterNew(m_vClusters.size(), cone.iMaxLevel);
			clusterNew.vConeIds.push_back(cone.iId);
			computeWorkLoad(clusterNew);
			//if(clusterNew.nNodes < 3) //prevent to cut small cone;
				//continue;
			if(workloadNewMax < clusterNew.iWorkload)
				workloadNewMax = clusterNew.iWorkload;
			m_vClusters.push_back(move(clusterNew));
			res++;

			//delete the cone from the cluster
			m_vClusters[clusterId].vConeIds[i] = -1;
		}
	}

	if(res) //PiCut happened!
	{
		ylog("In partClusterByPiCut(): %d PiCut happened for cluster[%d]\n", res, clusterId);
		int maxLev = 0;
		sort(m_vClusters[clusterId].vConeIds.begin(), m_vClusters[clusterId].vConeIds.end(), [](int32_t lhs, int32_t rhs){ return lhs > rhs; });
		for(int i = 0; i < m_vClusters[clusterId].vConeIds.size(); i++)
		{
			if(m_vClusters[clusterId].vConeIds[i] == -1) //delete cut cones
			{
				m_vClusters[clusterId].vConeIds.resize(i);
				break;
			}
			auto& cone = m_vCones[m_vClusters[clusterId].vConeIds[i]];
			if(maxLev < cone.iMaxLevel) //update maxlevel
				maxLev = cone.iMaxLevel;
		}
		if(m_vClusters[clusterId].iMaxLevel > maxLev) //critical cone is parted away
			m_vClusters[clusterId].iMaxLevel = maxLev;
		computeWorkLoad(m_vClusters[clusterId]);
		m_iMaxClusterWorkLoad = m_vClusters[clusterId].iWorkload; //assume the current cluster is still the biggest one.
		if(m_iMaxClusterWorkLoad <  workloadNewMax) //cut-out cone may be the biggest cluster?
			m_iMaxClusterWorkLoad = workloadNewMax;
		if(m_iMaxClusterWorkLoad < m_vClusters[clusterId + 1].iWorkload) //After PI-Cut, next cluster may be bigger?
			m_iMaxClusterWorkLoad = m_vClusters[clusterId + 1].iWorkload; 
	}

	if(!res)
		ylog("In partClusterByPiCut(): No PiCut happened for cluster[%d]\n", clusterId);
	return res;
}

void MetisAig::visitConeForEdgeWight(int32_t nodeId, int32_t coneId)
{
	Node& node = m_vNodes[nodeId];
	if(node.isPo())
	{
		visitConeForEdgeWight(node.getFanin0Id(), coneId);
		return ;
	}
	if(node.iIter == m_iGlobalIter)
		return ;
	if(node.isPi())
		return ;
	node.iIter = m_iGlobalIter;
	Cone& coneNode = m_vCones[node.iConeId];

	Node& fanin0 = m_vNodes[node.getFanin0Id()];
	Node& fanin1 = m_vNodes[node.getFanin1Id()];
	yassert((fanin0.iData != -1) && (fanin1.iData != -1) && (node.iData != -1));

	if(coneNode.iMaxLevel < m_iMaxLevel * CRITICAL_PATH_FACTOR / 100) //not critical
	{
		if(fanin0.isPi() && fanin0.nFanouts > 1)
			m_pMG->setEdgeWeight(fanin0.iData, node.iData, 1);
		else
			m_pMG->setEdgeWeight(fanin0.iData, node.iData, 100);

		if(fanin1.isPi() && fanin1.nFanouts > 1)
			m_pMG->setEdgeWeight(fanin1.iData, node.iData, 1);
		else
			m_pMG->setEdgeWeight(fanin1.iData, node.iData, 100);
	}
	else //critical path 
	{
		if(fanin0.isPi() && fanin0.nFanouts > 1)
			m_pMG->setEdgeWeight(fanin0.iData, node.iData, 1);
		else if(m_vCones[fanin0.iConeId].iMaxLevel > coneNode.iMaxLevel)
			m_pMG->setEdgeWeight(fanin0.iData, node.iData, 100);

		if(fanin1.isPi() && fanin1.nFanouts > 1)
			m_pMG->setEdgeWeight(fanin1.iData, node.iData, 1);
		else if(m_vCones[fanin1.iConeId].iMaxLevel > coneNode.iMaxLevel)
			m_pMG->setEdgeWeight(fanin1.iData, node.iData, 100);

	}

	visitConeForEdgeWight(fanin0.id, coneId);
	visitConeForEdgeWight(fanin1.id, coneId);

	return;


}

//划分AIG
//调用 decideNumParts()函数，决定 AIG 图的分区数量。如果返回值为 -1，则使用 Metis 算法进行分区，否则使用确定好的分区数量。
int32_t MetisAig::partitionAig()
{
	// int32_t nParts = decideNumParts();
	int32_t nParts = m_vClusters.size();
// 	if(nParts == -1)
// 	{
// 		//Metis is used. Prepare weight vectors
// 		m_pMG->initWeightVector();
// 		setNextIter();
// 		for(auto& cone : m_vCones)
// 		{
// 			visitConeForEdgeWight(cone.iPoId, cone.iId);
// 			//对 AIG 图中的节点进行遍历，并计算节点间的边权重和节点权重
// #if 0
// 			for(auto& edge : cone.vBoundaryEdges)
// 			{
				
// 				Node fanin = m_vNodes[edge.iFaninId];
// 				Node fanout = m_vNodes[edge.iFanoutId];
// 				if(fanout.isPo())
// 					continue;
// 				if(fanin.isPi())
// 					m_pMG->setEdgeWeight(fanin.iData, fanout.iData, 1);
// 				else 
// 					m_pMG->setEdgeWeight(fanin.iData, fanout.iData, 10);
// 			}
// #endif
// 		}
// 		for(auto& node : m_vNodes)
// 		{
// 			if(node.isPi() || node.isPo())
// 			{
// 				if(node.iData != -1)
// 					m_pMG->setNodeWeight(node.iData, 1);
// 			}
// 			else
// 				if(node.iData != -1)
// 					m_pMG->setNodeWeight(node.iData, 100);
// 		}
// 		//调用 partGraphByMetis() 方法进行分区
// 		m_pMG->partGraphByMetis(METIS_N_PART);
// 		return METIS_N_PART;
// 	}
	//yassert(nParts * m_iMaxClusterWorkLoad <= m_iTotalWorkLoad);

	//如果不适用metis，直接将 AIG 图分为指定数量的分区。
	vector<Partition> partitions(nParts);
	for(auto& cluster : m_vClusters)
	{
		partitions[0].addCluster(cluster);
		sort(partitions.begin(), partitions.end(), [](const Partition& lhs, const Partition& rhs){return lhs.iWorkload < rhs.iWorkload;});
	}

	/*
	for(auto& part : partitions)
	{
		ylog("This part: workload = %d\tnNodes = %d\n", part.iWorkload, part.nNodes);
	}
	*/

	for(int i = 0; i < partitions.size(); i++)
	{
		auto& part = partitions[i];
		for(auto clusterId : part.vClusterIds)
		{
			auto& cluster = m_vClusters[clusterId];
			setGraphPartition(cluster, i);
		}
	}

	m_pMG->setPoPart();
	m_pMG->setNumParts(nParts);
	return nParts;
}

void MetisAig::setGraphPartition(Cluster& cluster, int32_t partId)
{
	setNextIter();
	for(auto coneId : cluster.vConeIds)
		setNodePartition(m_vCones[coneId].iPoId, partId);
}

void MetisAig::setNodePartition(int32_t nodeId, int32_t partId)
{
	Node& node = m_vNodes[nodeId];
	if(node.isPo())
	{
		setNodePartition(node.getFanin0Id(), partId);
		return;
	}
	if(node.isPi())
	{
		if(node.iData != -1)
			m_pMG->setNodePart(node.iData, partId);
		return;
	}
	if(node.iIter == m_iGlobalIter)
		return; 
	node.iIter = m_iGlobalIter;
	if(node.iData != -1)
		m_pMG->setNodePart(node.iData, partId);
	setNodePartition(node.getFanin0Id(), partId);
	setNodePartition(node.getFanin1Id(), partId);
	return; 
}

/**Function*************************************************************

  Synopsis    [Merge several mapped networks.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
Abc_Ntk_t *Abc_NtkMerge(Abc_Ntk_t *pNtk, Vec_Ptr_t *pSubNtksOld, Vec_Ptr_t *pSubNtksNew)
{
	Abc_Ntk_t *pNtkRes, *pSubNtk;
	Abc_Obj_t *pObj, *pObjNew, *pFanin, *pConst1, *pConst0;
	Nm_Man_t *pManName;
	int iObjId;
	int i, j, k;
	pConst0 = pConst1 = NULL;

	yassert(Vec_PtrSize(pSubNtksNew) == Vec_PtrSize(pSubNtksOld));
	if (Vec_PtrSize(pSubNtksNew) == 0)
		return NULL;
	pNtkRes = Abc_NtkStartFrom(pNtk, ABC_NTK_LOGIC, ABC_FUNC_AIG);
	
	
	pManName = pNtkRes->pManName;
	// Clone information about newly generated PI/POs during partition process from pSubNtksOld to pSubNtksNew
	Vec_PtrForEachEntry(Abc_Ntk_t *, pSubNtksOld, pSubNtk, i)
    //for ( i = 0; (i < Vec_PtrSize(pSubNtksOld)) && (((pSubNtk) = (Type)Vec_PtrEntry(pSubNtksOld, i)), 1); i++ )
	{
		Abc_NtkForEachPi(pSubNtk, pObj, j)
    	//for ( i = 0; (i < Abc_NtkPiNum(pSubNtk)) && (((pObj) = Abc_NtkPi(pSubNtk, i)), 1); i++ )
		{
			if (pObj->fMarkA == 1)
			{
				pObj->pCopy->fMarkA = 1;
				pObj->pCopy->pData = Abc_ObjCopy(static_cast<Abc_Obj_t*>(pObj->pData));
				yassert(pObj->pCopy->pData != NULL);
			}
		}
		Abc_NtkForEachPo(pSubNtk, pObj, j)
		{
			if (pObj->fMarkA == 1)
				pObj->pCopy->fMarkA = 1;
		}
	}
	// Find exsiting CI/CO and Create new internal nodes in pNtkRes
	//Vec_PtrForEachEntry(Abc_Ntk_t *, pSubNtksNew, pSubNtk, i)
	Vec_PtrForEachEntry(Abc_Ntk_t *, pSubNtksNew, pSubNtk, i)
	{
		// Find new CI/COs in pNtkRes for corresponding CI/COs in pSubNtk
		Abc_NtkForEachCi(pSubNtk, pObj, j)
		{
			if (pObj->fMarkA != 1)
			{
				iObjId = Nm_ManFindIdByName(pManName, Abc_ObjName(pObj), Abc_ObjType(pObj) );
				if(iObjId == -1)
				{
					iObjId = Nm_ManFindIdByName(pManName, Abc_ObjName(pObj),  ABC_OBJ_BO );
				}
				pObj->pCopy = Abc_NtkObj(pNtkRes, iObjId);
			}
		}
		Abc_NtkForEachCo(pSubNtk, pObj, j)
		{
			if (pObj->fMarkA != 1)
			{
				iObjId = Nm_ManFindIdByName(pManName, Abc_ObjName(pObj), Abc_ObjType(pObj) );
				if(iObjId == -1)
				{
					iObjId = Nm_ManFindIdByName(pManName, Abc_ObjName(pObj),  ABC_OBJ_BI );
				}
				pObj->pCopy = Abc_NtkObj(pNtkRes, iObjId);
			}
		}
		// Create internel logic nodes
		Abc_NtkForEachNode(pSubNtk, pObj, j)
		{
			yassert(pObj->pData != NULL);
			if (Abc_ObjFaninNum(pObj) == 0) //ymc: take complement info into consideration
			{
				if (Hop_IsComplement(static_cast<Hop_Obj_t*>(pObj->pData)))
					if(pConst0)
						pObj->pCopy = pConst0;
					else
						pObj->pCopy = pConst0 = Abc_NtkCreateNodeConst0(pNtkRes);
				else
					if(pConst1)
						pObj->pCopy = pConst1;
					else
						pObj->pCopy = pConst1 = Abc_NtkCreateNodeConst1(pNtkRes);
				continue;
			}
			pObjNew = Abc_NtkCreateNode(pNtkRes);
			pObjNew->pData = Hop_Transfer(static_cast<Hop_Man_t*>(pSubNtk->pManFunc), static_cast<Hop_Man_t*>(pNtkRes->pManFunc), static_cast<Hop_Obj_t*>(pObj->pData), Abc_ObjFaninNum(pObj));
			pObj->pCopy = pObjNew;
		}
	}
	//Finish Connections
	Vec_PtrForEachEntry(Abc_Ntk_t *, pSubNtksNew, pSubNtk, i)
	{
		Abc_NtkForEachCo(pSubNtk, pObj, j)
		{
			pFanin = Abc_ObjFanin0(pObj);
			if (pObj->fMarkA != 1) //cut-caused PI/POs don't need connection
			{
				if (pFanin->fMarkA != 1)
					Abc_ObjAddFanin(pObj->pCopy, Abc_ObjNotCond(pFanin->pCopy, pObj->fCompl0));
				else //if CO's fanin is a cut-caused CI
				{
					if(Abc_ObjIsPo(static_cast<Abc_Obj_t*>(pFanin->pData)))
						Abc_ObjAddFanin(pObj->pCopy, Abc_ObjNotCond(Abc_ObjCopy(Abc_ObjFanin0(static_cast<Abc_Obj_t*>(pFanin->pData))), Abc_ObjFaninC0(static_cast<Abc_Obj_t*>(pFanin->pData)) ^ Abc_ObjFaninC0(pObj)));
					else if(Abc_ObjIsPi(static_cast<Abc_Obj_t*>(pFanin->pData)))
						Abc_ObjAddFanin(pObj->pCopy, Abc_ObjNotCond(Abc_ObjCopy(static_cast<Abc_Obj_t*>(pFanin->pData)), Abc_ObjFaninC0(pObj)));
					else
						yassert(0);
				}
			}
		}
		Abc_NtkForEachNode(pSubNtk, pObj, j)
		{
			Abc_ObjForEachFanin(pObj, pFanin, k)
			{
				if (pFanin->fMarkA != 1)
					Abc_ObjAddFanin(pObj->pCopy, pFanin->pCopy);
				else
				{
					if(Abc_ObjIsPo(static_cast<Abc_Obj_t*>(pFanin->pData)))
					{
						if (Abc_ObjFaninC0(static_cast<Abc_Obj_t*>(pFanin->pData)) == 1)
						{
							pObj->pData = Hop_Complement(static_cast<Hop_Man_t*>(pNtkRes->pManFunc), static_cast<Hop_Obj_t*>(pObj->pData), k);
						}
						Abc_ObjAddFanin(pObj->pCopy, Abc_ObjCopy(Abc_ObjFanin0(static_cast<Abc_Obj_t*>(pFanin->pData))));
					}
					else if(Abc_ObjIsPi(static_cast<Abc_Obj_t*>(pFanin->pData)))
						Abc_ObjAddFanin(pObj->pCopy, Abc_ObjNotCond(Abc_ObjCopy(static_cast<Abc_Obj_t*>(pFanin->pData)), Abc_ObjFaninC0(pObj)));
					else
						yassert(0);
				}
			}
		}
	}

	//zkx
	// for ( i = 0; (i < Vec_PtrSize(pSubNtksOld)) && (((pNtk) = (Abc_Ntk_t *)Vec_PtrEntry(pSubNtksOld, i)), 1); i++ ){
	// //for(auto pNtk : pSubNtksNew){

	// 	char filename[256];
	// 	const char* outputFolder ="/home/kxzhu/partition/abc/test";
	// 	snprintf(filename, sizeof(filename),"%s/Mappednetwork_%d.blif",outputFolder, i);
	// 	printf("Writing file: %s\n", filename);

	// 	Abc_Ntk_t* pNtkNew = Abc_NtkToNetlist(pNtk);
	// 	//cout<<"Now in iteration!"<<endl;

	// 	Io_WriteBlif(pNtkNew,filename,1,0,0);
		
	// }
	
	// wxx: Network check. Attention: some of PI->pData are not empty!
	// if ( !Abc_NtkCheck( pNtkRes ) )
	// {
	// 	printf( "Abc_NtkPartition: The network check has failed.\n" );
	// 	Abc_NtkCleanMarkA(pNtkRes);
	// 	Abc_NtkDelete( pNtkRes );
	// 	return NULL;
	// }
	//	Abc_NtkShow( pSubNtksNew->pArray[0], 0, 0, 1);

	return pNtkRes;
}

static Abc_Obj_t* Abc_ObjExtractSubNtk_rec(Abc_Obj_t* pObj, Abc_Ntk_t* pNtkNew)
{
	if(Abc_NodeIsTravIdCurrent(pObj))
		return pObj->pCopy;
	Abc_NodeSetTravIdCurrent(pObj);
	if(Abc_ObjIsCi(pObj))
	{
		pObj->pCopy = Abc_NtkCreatePi(pNtkNew);
		Abc_ObjAssignName(pObj->pCopy, Abc_ObjName(pObj), NULL);
		return pObj->pCopy;
	}
	if(Abc_ObjIsNode(pObj))
	{
		Abc_Obj_t* pFanin0New = Abc_ObjExtractSubNtk_rec(Abc_ObjFanin0(pObj), pNtkNew);
		Abc_Obj_t* pFanin1New = Abc_ObjExtractSubNtk_rec(Abc_ObjFanin1(pObj), pNtkNew);
		pObj->pCopy = Abc_AigAnd(static_cast<Abc_Aig_t*>(pNtkNew->pManFunc), Abc_ObjNotCond(pFanin0New, Abc_ObjFaninC0(pObj)), Abc_ObjNotCond(pFanin1New, Abc_ObjFaninC1(pObj)));
		return pObj->pCopy;
	}
	if(Abc_ObjIsCo(pObj))
	{
		Abc_Obj_t* pFaninNew = Abc_ObjExtractSubNtk_rec(Abc_ObjFanin0(pObj), pNtkNew);
		pObj->pCopy = Abc_NtkCreatePo(pNtkNew);
		Abc_ObjAssignName(pObj->pCopy, Abc_ObjName(pObj), NULL);
		Abc_ObjAddFanin(pObj->pCopy, pFaninNew);
		return pObj->pCopy;
	}
	yassert(0);
}


Abc_Ntk_t* Abc_NtkExtractCriticalPath(Abc_Ntk_t* pNtk, Abc_Obj_t* pCo)
{
	yassert(pNtk->ntkType == ABC_NTK_STRASH);	
	yassert(pCo->pNtk == pNtk);
	yassert(Abc_ObjIsCo(pCo));
	Abc_Ntk_t* pNtkRes = Abc_NtkAlloc(pNtk->ntkType, pNtk->ntkFunc, 1);
	pNtkRes->nConstrs = pNtk->nConstrs;
	pNtkRes->nBarBufs = pNtk->nBarBufs;
	// duplicate the name and the spec
	pNtkRes->pName = Extra_UtilStrsav(pNtk->pName);
	pNtkRes->pSpec = Extra_UtilStrsav(pNtk->pSpec);

	Abc_NtkIncrementTravId(pNtk);
	Abc_ObjExtractSubNtk_rec(pCo, pNtkRes);
	printf("Num of Pi: %d\n", Abc_NtkPiNum(pNtkRes));
	printf("Num of Po: %d\n", Abc_NtkPoNum(pNtkRes));
	printf("Num of AND Node: %d\n", Abc_NtkNodeNum(pNtkRes));
	return pNtkRes;
}

Abc_Obj_t* Abc_NtkPickCriticalPo(Abc_Ntk_t* pNtk)
{
	printf("Now pick Po\n");
	yassert(pNtk->ntkType == ABC_NTK_LOGIC);
	int levelMax = Abc_NtkLevel(pNtk);
	printf("levelMax = %d\n", levelMax);
	Abc_Obj_t* pObj;
	int i;
	Abc_NtkForEachCo(pNtk, pObj, i)
	{
		Abc_Obj_t* pDriver = Abc_ObjFanin0(pObj);
		if(pDriver->Level == levelMax)
		{
			printf("Name of Po: %s\n", Abc_ObjName(Abc_ObjFanout0(pObj)));
			printf("Type of Co's fanout: %d\n", Abc_ObjType(Abc_ObjFanout0(pObj)));
			printf("Level of Po: %d\n", levelMax);
			printf("i of Po : %d\n", i);
			return pObj;
		}
	}
	return NULL;
}

} //for namespace

