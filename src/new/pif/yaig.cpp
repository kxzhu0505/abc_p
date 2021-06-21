#include "yaig.h"

namespace ymc {

int32_t Yaig::addAndNode(Literal fanin0, Literal fanin1)
{
	yassert(fanin0.getId() >= 2); //0 is for NONE; 1 is for CONST
	yassert(fanin0.getId() < m_vNodes.size());
	yassert(fanin1.getId() < m_vNodes.size());
	if (fanin1.getId()) //sort
		if (fanin1.getId() < fanin0.getId())
			swap(fanin0, fanin1);
	NodeAttribute tmpNodeAttr(fanin0, fanin1);
	auto iter = m_strashTable.find(tmpNodeAttr.toUint64());
	if (iter == m_strashTable.end())
	{
		int32_t newId = m_vNodes.size();
		Node node(tmpNodeAttr, newId);
		m_vNodes.push_back(move(node));
		m_strashTable.emplace(tmpNodeAttr.toUint64(), newId);
		m_vNodes[fanin0.getId()].nFanouts++;
		m_nEdges++;
		if (fanin1.getId() >= 2)
		{
			m_vNodes[fanin1.getId()].nFanouts++;
			m_nEdges++;
		}
		return newId;
	}
	else
	{
		return m_vNodes[iter->second].id;
	}
}

int32_t Yaig::addAndNode(int32_t fanin0Id, bool isC0, int32_t fanin1Id, bool isC1)
{
	return addAndNode(Literal(fanin0Id, isC0), Literal(fanin1Id, isC1));
}

int32_t Yaig::addPi()
{
	int32_t newId = m_vNodes.size();
	Node node(NodeAttribute(0), newId);
	m_vNodes.push_back(move(node));
	m_vPis.push_back(newId);
	return newId;
}

int32_t Yaig::addPo(const Literal& fanin)
{
	int32_t curN = m_vNodes.size();
	int32_t newId = addAndNode(fanin, 0);
	if(newId == curN)
		m_vPos.push_back(newId);
	return newId;
}

int32_t Yaig::addPo(const int32_t faninId, const bool isC)
{
	return addPo(Literal(faninId, isC));
}

void Yaig::addSubGraph(const Graph& graph) //For Function Library
{
	graph.check();
	int32_t nNodes = graph.getNodeNum();
	yassert(nNodes > 0);
	int32_t nPis = graph.getPiNum();
	vector<int32_t> vNodeMap(nNodes); //NodeMap: nodes of subGraph (graph) -> nodeIds in Yaig
	int32_t cnt = 0, fanin0g, fanin1g;
	bool isC0g, isC1g;

	for (int32_t i = m_vPis.size(); i < nPis; i++) //make sure that there are enough PIs in AIG
		addPi();

	for (int32_t i = 0; i < nNodes; i++)
	{
		if (graph.isNodePi(i))
			vNodeMap[i] = m_vPis[cnt++];
	}	
	yassert(cnt == nPis);

	for(int32_t i = 0; i < nNodes ; i++ )
	{
		if (graph.isNodePo(i))
		{
			fanin0g = graph.getNodeFanin0Id(i);
			isC0g = graph.isNodeFanin0C(i);
			vNodeMap[i] = addPo(vNodeMap[fanin0g], isC0g);
		}
		else if (graph.isNodeAND(i)) //for AND
		{
			fanin0g = graph.getNodeFanin0Id(i);
			isC0g = graph.isNodeFanin0C(i);
			fanin1g = graph.getNodeFanin1Id(i);
			isC1g = graph.isNodeFanin1C(i);
			vNodeMap[i] = addAndNode(vNodeMap[fanin0g], isC0g,
				vNodeMap[fanin1g], isC1g);
		}
	}

	for(int i = 0; i < vNodeMap.size(); i++)
		m_vNodes[vNodeMap[i]].iData = i;
	return;
}

shared_ptr<Graph> Yaig::toGraph()
{
	yassert(!m_vNodes.empty());
	//TODO: const0?
	int32_t nNodes = m_vNodes.size() - ID_OFFSET;
	shared_ptr<Graph> spg = make_shared<Graph>(nNodes, m_nEdges, m_vNodes);
	return spg;
}

void Yaig::clear()
{
	vector<Node>().swap(m_vNodes);
	unordered_map<uint64_t, int32_t>().swap(m_strashTable);
	vector<int32_t>().swap(m_vPis);
	vector<int32_t>().swap(m_vPos);
	m_nEdges = 0;
	m_iGlobalIter = 0;
	m_iMaxLevel = 0;
}

void Yaig::check()
{
	ylog("Now in yaig::check()\n");
	ylog("total number of nodes = %ld\n", m_vNodes.size() - 2);
	ylog("total number of edges = %d\n", m_nEdges);
	ylog("total number of PI = %ld\n", m_vPis.size());
	ylog("total number of PO = %ld\n", m_vPos.size());
	int32_t cnt = 0;
	int32_t i;
	for (i = 2; i < m_vNodes.size(); i++)
		cnt += m_vNodes[i].nFanouts;
	ylog("total number of fanouts = %d\n", cnt);
	for (i = 2, cnt = 0; i < m_vNodes.size(); i++)
	{
		if (m_vNodes[i].isAnd())
			cnt += 2;
		else if (m_vNodes[i].isPo())
			cnt++;
	}
	ylog("total number of fanins = %d\n", cnt);
	for (auto node : m_vNodes)
	{
		yassert(node.getFanin0Id() == 0 || node.id > node.getFanin0Id()); //DFS
		yassert(node.getFanin0Id() == 0 || node.id > node.getFanin1Id());
	}
#if 0
	ylog("print all nodes below:\n");
	for (auto i : m_vNodes)
		printOneNode(i);	
#endif

}

void Yaig::printOneNode(Node& node)
{
	ylog("node id: %d\tnFanouts: %d\n", node.id-ID_OFFSET, node.nFanouts);
	ylog("fanin0id: %d, isC0: %d\n", node.getFanin0Id(), node.isC0());
	ylog("fanin1id: %d, isC1: %d\n", node.getFanin1Id(), node.isC1());
	ylog("iLevel: %d\tiData: %d\tiIter: %d\tiConeId: %d\tiNCUTs = %d\n", node.iLevel, node.iData, node.iIter, node.iConeId, node.iNCuts);
}

int32_t Yaig::computeAllLevel()
{
	ylog("Now in computeAllLevel()\n");
	int level;
	setNextIter();
	for(auto po : m_vPos)
	{
		level = computeNodeLevel(m_vNodes[po].getFanin0Id());
		m_vNodes[po].iLevel = level;
		if(level > m_iMaxLevel)
			m_iMaxLevel = level;
	}
	ylog("MaxLevel: %d\n", m_iMaxLevel);
	return m_iMaxLevel;
}

int32_t Yaig::computeNodeLevel(int32_t id)
{
	int32_t res = 1;
	int32_t level0, level1;
	yassert(!m_vNodes[id].isPo());
	if(m_vNodes[id].isPi())
		return 0;
	if(m_vNodes[id].iIter == getCurrentIter())
		return m_vNodes[id].iLevel;
	m_vNodes[id].iIter = getCurrentIter();
	level0 = computeNodeLevel(m_vNodes[id].getFanin0Id());
	level1 = computeNodeLevel(m_vNodes[id].getFanin1Id());
	res += (level0 < level1) ? level1 : level0;
	m_vNodes[id].iLevel = res;
	return res;
}



void Graph::print()
{
	ylog("nNodes = %d\n", m_nNodes);
	ylog("vNodeIndices: ");
	for (auto i : m_vNodeIndices)
		ylog("%d ", i);
	ylog("\n");
	ylog("vNodeIsPi: ");
	for (auto i : m_vNodeIsPi)
		ylog(i ? "1 " : "0 ");
	ylog("\n");
	ylog("nEdges = %ld\n", m_vEdges.size());
	ylog("vEdges: ");
	for (auto i : m_vEdges)
		ylog("%d ", i);
	ylog("\n");
	ylog("vEdgeIsC: ");
	for (auto i : m_vEdgeIsC)
		ylog(i ? "1 " : "0 ");
	ylog("\n");
}

void Graph::check() const
{
	yassert(m_nNodes == m_vNodeIndices.size() - 1);
	yassert(m_nNodes == m_vNodeIsPi.size());
	yassert(m_vEdges.size() == m_vEdgeIsC.size());
	for (int32_t iNode = 0; iNode < m_nNodes; iNode++)
	{
		if (m_vNodeIsPi[iNode]) //iNode is PI
			for (int32_t j = m_vNodeIndices[iNode]; j < m_vNodeIndices[iNode + 1]; j++)
				yassert(iNode < m_vEdges[j]);
		else if (m_vNodeIndices[iNode + 1] == m_vNodeIndices[iNode] + 1) //iNode is PO
			yassert(iNode > m_vEdges[m_vNodeIndices[iNode]]);
		else //iNode is AND Node
		{
			int32_t nNeighbor = m_vNodeIndices[iNode + 1] - m_vNodeIndices[iNode];
			yassert(nNeighbor > 2);
			yassert(iNode > m_vEdges[m_vNodeIndices[iNode]]);
			yassert(iNode > m_vEdges[m_vNodeIndices[iNode] + 1]);
			int32_t iFanout = m_vNodeIndices[iNode] + 2;
			for (int32_t iEdge = iFanout; iEdge < m_vNodeIndices[iNode + 1]; iEdge++)
				yassert(iNode < m_vEdges[iEdge]);
		}
	}
}

void Graph::init(int32_t nNodes, int32_t nEdges, int32_t nParts)
{
	m_nNodes = nNodes;
	m_nCon = 1;
	m_nPis = 0;
	m_vNodeIndices.resize(nNodes + 1);
	m_vEdges.resize(nEdges * 2);
	m_vEdgeIsC.resize(nEdges * 2);
	m_vNodeIsPi.resize(nNodes);
	m_nParts = nParts;
	m_nCutEdges = 0;
	m_vPartition.resize(nNodes);
	m_vNodeIndices[0] = 0;
}

Graph::Graph(const int32_t nNodes, const int32_t nEdges, const vector<Node>& vNodes)
{
	int32_t fanin0, fanin1, edge;
	vector<int32_t> iEdge(nNodes, 0); //iter for edges of each node
	int32_t offset = vNodes.size() - nNodes;
	yassert(offset >= 0);

	init(nNodes, nEdges, 0); //resize internal vectors

	for (int32_t iNode = 0, i = offset; i < vNodes.size(); i++)
	{
		if (vNodes[i].isPi())
		{
			m_vNodeIndices[iNode + 1] = m_vNodeIndices[iNode] + vNodes[i].nFanouts;
			m_vNodeIsPi[iNode] = true;
			m_nPis++;
		}
		else if (vNodes[i].isPo())
		{
			m_vNodeIndices[iNode + 1] = m_vNodeIndices[iNode] + 1;
			m_vNodeIsPi[iNode] = false;
			fanin0 = vNodes[i].getFanin0Id() - offset; //NodeId in graph

			edge = m_vNodeIndices[iNode] + iEdge[iNode]++; //fanin for current PO node
			m_vEdges[edge] = fanin0;
			m_vEdgeIsC[edge] = vNodes[i].isC0();

			edge = m_vNodeIndices[fanin0] + iEdge[fanin0]++; //corresponding fanout 
			m_vEdges[edge] = iNode;
			m_vEdgeIsC[edge] = vNodes[i].isC0();
		}
		else
		{
			m_vNodeIndices[iNode + 1] = m_vNodeIndices[iNode] + 2 + vNodes[i].nFanouts;
			m_vNodeIsPi[iNode] = false;
			fanin0 = vNodes[i].getFanin0Id() - offset;
			fanin1 = vNodes[i].getFanin1Id() - offset;

			//for fanin0
			edge = m_vNodeIndices[iNode] + iEdge[iNode]++; //fanin0 for current AND node
			m_vEdges[edge] = fanin0;
			m_vEdgeIsC[edge] = vNodes[i].isC0();
			edge = m_vNodeIndices[fanin0] + iEdge[fanin0]++;
			m_vEdges[edge] = iNode;
			m_vEdgeIsC[edge] = vNodes[i].isC0();

			//for fanin1
			edge = m_vNodeIndices[iNode] + iEdge[iNode]++; //fanin1 for current AND node
			m_vEdges[edge] = fanin1;
			m_vEdgeIsC[edge] = vNodes[i].isC1();
			edge = m_vNodeIndices[fanin1] + iEdge[fanin1]++;
			m_vEdges[edge] = iNode;
			m_vEdgeIsC[edge] = vNodes[i].isC1();
		}
		iNode++;
	}

}

} //for namespaces