#pragma once

extern "C" {
#include "base/abc/abc.h"
#include "map/if/if.h"
#include "map/mio/mio.h"
#include "base/main/main.h"
}
#include<sys/time.h>
#include<metis.h>
#include "yaig.h"
#include <map>
#include <utility>

namespace ymc {

class MetisGraph:public Graph //for partitioning
{
public:
    MetisGraph() = delete;
    ~MetisGraph() = default;
    MetisGraph(Abc_Ntk_t* pNtk, bool fMetis);
    void metisGraphInit(int32_t nNodes, bool fMetis);
    int createGraphFromNtk(Abc_Ntk_t* pNtk, bool fMetis);
    void initWeightVector(){
        yassert(m_nNodes);
        yassert(!m_vEdges.empty());
        m_vNodeWeights.resize(m_nNodes, 1);
        m_vEdgeWeights.resize(m_vEdges.size(), WALL_EDGE_WEIGHT);
    }

    int partGraphByMetis(int32_t nParts);
    Abc_Ntk_t* initOneSubNtk(int index);
    int createSubNtksFromPartition(vector<Abc_Ntk_t*>& vSubNtks);

    void tmp();
    int computeNodeWeightByLevel(int leftChildLevel, int rightChildLevel);
    int computeNodeWeightByDegree(int degree);
    int computeEdgeWeight(int faninLevel, int faninReqTime, int level, int reqTime);

    void setEdgeWeight(int fanin, int fanout, int weight);
    void setNodeWeight(int nodeId, int weight)
    {
        yassert(nodeId < m_nNodes);
        yassert(nodeId < m_vNodeWeights.size());
        m_vNodeWeights[nodeId] = weight;
    }
    void setNodePart(int id, int ipart) {
        if(m_vPartition[id] == -1) {
            m_vPartition[id] = ipart; 
            return;
        }
        else {
            if(m_vPartition[id] == ipart)
                return;
            else
                if(isNodeAND(id)) {
                    ylog("Node[%d] has multi-partition: %d and %d\n", id, ipart, m_vPartition[id]);
                    yassert(0);
                }
                else
                    return;
        }
    }
    void setPoPart();

private:
    const int MAX_NODE_WEIGHT = 1; 
    const int MAX_EDGE_WEIGHT_FOR_NODE = 150; 
    const int MIN_EDGE_WEIGHT_FOR_NODE = 10; 
    const int WALL_EDGE_WEIGHT = 10000; //wall to prevent metis to cut the edge 
    const int WALL_THRESHOLD = 20; //threshold of slack to generate wall
    const int METIS_UFACTOR = 30;
    Abc_Ntk_t* m_pNtkOrigin;
    vector<Abc_Obj_t*> m_vpObjs; //nodeId -> pObj in origin Ntk
    vector<Abc_Obj_t*> m_vpObjsNew; //nodeId -> pObj in SubNtks
    vector<int32_t> m_vNodeWeights;
    vector<int32_t> m_vEdgeWeights;
};

class Edge {
public:
    Edge(int32_t in, int32_t out):iFaninId(in), iFanoutId(out){};
	int32_t iFaninId;
	int32_t iFanoutId;
};

class Cone {
public:
    Cone(int32_t lev, int32_t nodeId) : iMaxLevel(lev), iPoId(nodeId), iId(0){};
	int32_t iMaxLevel;
	int32_t iPoId;
    int32_t iId;
	vector<Edge> vBoundaryEdges;
    set<int32_t> sAdjacentConeId;
};

class Cluster {
public: 
    Cluster() = default;
    Cluster(int32_t id, int32_t lev):iId(id), iMaxLevel(lev), iWorkload(0), nNodes(0), iPartitionId(-1) {};
    int32_t iId;
	int32_t iMaxLevel;
	int32_t iWorkload;
    int32_t nNodes; //excluding PI/PO
    int32_t iPartitionId;
	vector<int32_t> vConeIds;
};

class Partition {
public:
    Partition(): iWorkload(0), nNodes(0) {};
    void addCluster(Cluster& cluster){
        iWorkload += cluster.iWorkload;
        nNodes += cluster.nNodes;
        vClusterIds.push_back(cluster.iId);
    };
    int32_t iWorkload;
    int32_t nNodes;
    vector<int32_t> vClusterIds;
};

class MetisAig : public Yaig //for analysing the original NTK
{
public:
    MetisAig() : m_pMG(NULL) {} ;
    ~MetisAig() = default;
    void bindGraph(MetisGraph* pmg);

    void parseAig();
	void visitAllFaninFromNode(int32_t nodeId, Cone& cone); //set node.iIter equal to m_iGlobalIter
    void findBoundaryEdges(Cone& cone);
    void findBoundaryEdges_rec(Cone& cone, int32_t nodeId, int32_t fCovered);
    void computeWorkLoad(Cluster& cluster);
    int32_t computeWorkLoad_rec(int32_t nodeId, int32_t& workload, int32_t& nNodes);
    int32_t decideNumParts();
    int32_t partBiggestClusterByPICut(int32_t clusterId);
    int32_t partBiggestCluster(int32_t clusterId, int32_t workLoadLimit);
    void cutOneBoundaryEdge(Edge& edge);
    int32_t tryPart(); //A-B-C
    int32_t tryPart2(); //from critical
    int32_t partitionAig();
    void setGraphPartition(Cluster& cluster, int32_t partId);
    void setNodePartition(int32_t nodeId, int32_t partId);
    void visitConeForEdgeWight(int32_t nodeId, int32_t coneId);

    void printCones();
    void printOneCone(int32_t coneId);
    void printOneCone_rec(int32_t nodeId);
    void printClusters();
    void checkClusters();
    void checkNodeCluster(int32_t nodeId, int32_t clusterId);

private:
    const int32_t PI_WORK_LOAD = 1;
    const int32_t MAX_N_CUT = 8;
    const int32_t MIN_N_PART = 2;
    const int32_t MAX_N_PART = 20;
    const int32_t METIS_N_PART = 4;
    const int32_t CRITICAL_PATH_FACTOR = 50;

    MetisGraph* m_pMG;
    vector<Cone> m_vCones;
    vector<Cluster> m_vClusters;
    int32_t m_iTotalWorkLoad;
    int32_t m_iMaxClusterWorkLoad;

    vector<int> m_vConeId2ClusterId;
};


Abc_Ntk_t * Abc_NtkMerge( Abc_Ntk_t * pNtk, Vec_Ptr_t * pSubNtksOld, Vec_Ptr_t * pSubNtksNew);
Abc_Ntk_t* Abc_NtkExtractCriticalPath(Abc_Ntk_t* pNtk, Abc_Obj_t* pCos);
Abc_Obj_t* Abc_NtkPickCriticalPo(Abc_Ntk_t* pNtk);

} //for namespace
