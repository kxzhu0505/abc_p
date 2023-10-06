#pragma once
#include <iostream>
#include <cstdio>
#include <cassert>
#include <unordered_map>
#include <vector>
#include <set>
#include <memory>
#include <mutex>
#include <algorithm>
#include <numeric>
#include <thread>

#define ASSERT_DEBUG
#ifdef ASSERT_DEBUG
#define yassert(x) assert(x)
#else
#define yassert(x) 
#endif

#define LOG_DEBUG
#ifdef LOG_DEBUG
#define ylog(x...) printf(x)
#else
#define ylog(x...)
#endif

using namespace std;

namespace ymc {

class Literal {
public:
//isC 表示 isComplemented
	Literal(int32_t id, bool isC) :nodeId(id), isC(isC) {};
	Literal(int32_t in = 0) { fromUint32(in); }
	bool isComplemented() { return isC; }
	int32_t getId() { return nodeId; }
	int32_t toUint32() { return static_cast<int32_t>(isC) | (static_cast<int32_t>(nodeId) << 1); }
	void fromUint32(int32_t in) { isC = in & 1;	nodeId = in >> 1; }
	operator int32_t() { return toUint32(); }

	int32_t getId() const { return nodeId; }
	bool isNone() const { return nodeId == 0; }
	bool isConst0() const { return (nodeId == 1) && (isC == 1); }
	bool isConst1() const { return (nodeId == 1) && (isC == 0); }

private:
	unsigned nodeId : 31;
	unsigned isC : 1;
};

class NodeAttribute {
public:
	NodeAttribute(uint64_t in) { fromUint64(in); }
	NodeAttribute(Literal f0, Literal f1) : fanin0(f0), fanin1(f1) {};
	inline uint64_t toUint64() {
		return static_cast<uint64_t>(fanin0.toUint32()) | (static_cast<uint64_t>(fanin1.toUint32()) << 32);
	}
	inline void fromUint64(uint64_t in) {
		fanin0.fromUint32(static_cast<int32_t>(in));
		fanin1.fromUint32(static_cast<int32_t>(in >> 32));
	}
	operator uint64_t() { return toUint64(); }
	Literal getFanin0() const { return fanin0; }
	Literal getFanin1() const { return fanin1; }

private:
	Literal fanin0;
	Literal fanin1;
};

class Node {
public:
	Node(NodeAttribute attr, int32_t id) :attr(attr), id(id), nFanouts(0), nVisits(0), iLevel(-1), iNCuts(0), iData(-1), iIter(0), iConeId(-1) {};
	int32_t getFanin0Id() const { return attr.getFanin0().getId(); }
	int32_t getFanin1Id() const { return attr.getFanin1().getId(); }
	bool isC0() const { return attr.getFanin0().isComplemented(); }
	bool isC1() const { return attr.getFanin1().isComplemented(); }
	bool isPi() const { return attr.getFanin0().getId() == 0 && attr.getFanin1().getId() == 0; }
	bool isPo() const { return attr.getFanin0().getId() != 0 && attr.getFanin1().getId() == 0; }
	bool isAnd() const { return attr.getFanin0().getId() && attr.getFanin1().getId(); }
public:
	NodeAttribute attr;
	int32_t id;
	int32_t nFanouts; 
	int32_t nVisits; 
	int32_t iLevel; //logic height
	int32_t iNCuts; //for computation of workload
	int32_t iData; //used for corresponding node in MetisGraph
	uint32_t iIter; //used for iteration
	int32_t iConeId; //longest cone which this node belongs to
};

class Graph {
public:
	Graph() = default;
	Graph(const int32_t nNodes, const int32_t nEdges, const vector<Node>& vNodes); //init from yaig
	void init(int32_t nNodes, int32_t nEdges, int32_t nParts);
	void print();
	void check() const; //make sure the nodeId is bigger than its faninId, and smaller than its fanoutId
	void setNumParts(int32_t nParts) { m_nParts = nParts; }

	int32_t getNodeNum() const { return m_nNodes; }
	int32_t getPiNum() const { return m_nPis; }
	int32_t getNodeEdgeNum(int32_t nodeId) const { return m_vNodeIndices[nodeId + 1] - m_vNodeIndices[nodeId]; }
	bool isNodePi(int32_t nodeId) const { return m_vNodeIsPi[nodeId]; }
	bool isNodePo(int32_t nodeId) const { return !isNodePi(nodeId) && (getNodeEdgeNum(nodeId) == 1); }
	bool isNodeAND(int32_t nodeId) const { return !isNodePi(nodeId) && !isNodePo(nodeId); }
	int32_t getNodeFanin0Id(int32_t nodeId) const { yassert(!isNodePi(nodeId)); return m_vEdges[m_vNodeIndices[nodeId]]; }
	int32_t getNodeFanin1Id(int32_t nodeId) const { yassert(isNodeAND(nodeId)); return m_vEdges[m_vNodeIndices[nodeId] + 1]; }
	bool isNodeFanin0C(int32_t nodeId) const { yassert(!isNodePi(nodeId)); return m_vEdgeIsC[m_vNodeIndices[nodeId]]; }
	bool isNodeFanin1C(int32_t nodeId) const { yassert(isNodeAND(nodeId)); return m_vEdgeIsC[m_vNodeIndices[nodeId] + 1]; }

protected:
	int32_t m_nNodes;
	int32_t m_nCon;
	int32_t m_nPis;
	vector<int32_t> m_vNodeIndices; //size: N+1
	vector<int32_t> m_vEdges; //size: 2*M
	vector<bool> m_vEdgeIsC; //size: 2*M
	vector<bool> m_vNodeIsPi; //size: N
	int32_t m_nParts;
	int32_t m_nCutEdges;
	vector<int32_t> m_vPartition; //size: N
};

class Yaig {
public:
	const int32_t ID_OFFSET = 2;
	Yaig() : m_nEdges(0), m_iGlobalIter(0), m_iMaxLevel(0), m_tmp0(0), m_tmp1(0) 
		{m_vNodes.push_back(Node(0,0));m_vNodes.push_back(Node(0,1));}
	~Yaig() = default;
	int32_t addAndNode(Literal fanin0, Literal fanin1);
	int32_t addAndNode(int32_t fanin0Id, bool isC0, int32_t fanin1Id, bool isC1);
	int32_t addPi();
	int32_t addPo(const Literal& fanin);
	int32_t addPo(const int32_t faninId, const bool isC);
	void addSubGraph(const Graph& graph);
	shared_ptr<Graph> toGraph();
	void check();
	void clear();
	void printOneNode(Node& node);
	int32_t getNodeData(int32_t nodeId) const { return m_vNodes[nodeId].iData; }

	void setNextIter() { m_iGlobalIter++; }
	int32_t getCurrentIter() const { return m_iGlobalIter; }
	int32_t computeAllLevel();
	int32_t computeNodeLevel(int32_t id);

protected:
	unordered_map<uint64_t, int32_t> m_strashTable; //NodeAttribute -> NodeId
	vector<Node> m_vNodes; //v[0] means NONE, v[1] means CONST!!
	vector<int32_t> m_vPis;
	vector<int32_t> m_vPos;
	int32_t m_nEdges;

	uint32_t m_iGlobalIter;
	int32_t m_iMaxLevel;

	int32_t m_tmp0;
	int32_t m_tmp1;
};




}//for namespace