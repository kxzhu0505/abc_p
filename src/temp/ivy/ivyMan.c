/**CFile****************************************************************

  FileName    [ivyMan.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [And-Inverter Graph package.]

  Synopsis    [AIG manager.]

  Author      [Alan Mishchenko]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - May 11, 2006.]

  Revision    [$Id: ivy_.c,v 1.00 2006/05/11 00:00:00 alanmi Exp $]

***********************************************************************/

#include "ivy.h"

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    [Starts the AIG manager.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Ivy_Man_t * Ivy_ManStart()
{
    Ivy_Man_t * p;
    // start the manager
    p = ALLOC( Ivy_Man_t, 1 );
    memset( p, 0, sizeof(Ivy_Man_t) );
    // perform initializations
    p->Ghost.Id = -1;
    p->nTravIds = 1;
    p->fCatchExor = 1;
    // allocate arrays for nodes
    p->vPis = Vec_PtrAlloc( 100 );
    p->vPos = Vec_PtrAlloc( 100 );
    p->vBufs = Vec_PtrAlloc( 100 );
    p->vObjs = Vec_PtrAlloc( 100 );
    // prepare the internal memory manager
    Ivy_ManStartMemory( p );
    // create the constant node
    p->pConst1 = Ivy_ManFetchMemory( p );
//    p->pConst1->fPhase = 1;
    Vec_PtrPush( p->vObjs, p->pConst1 );
    p->nCreated = 1;
    // start the table
    p->nTableSize = 10007;
    p->pTable = ALLOC( int, p->nTableSize );
    memset( p->pTable, 0, sizeof(int) * p->nTableSize );
    return p;
}

/**Function*************************************************************

  Synopsis    [Duplicates the AIG manager.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Ivy_Man_t * Ivy_ManDup( Ivy_Man_t * p )
{
    Vec_Int_t * vNodes, * vLatches;
    Ivy_Man_t * pNew;
    Ivy_Obj_t * pObj;
    int i;
    // collect latches and nodes in the DFS order
    vNodes = Ivy_ManDfsSeq( p, &vLatches );
    // create the new manager
    pNew = Ivy_ManStart();
    // create the PIs
    Ivy_ManConst1(p)->pEquiv = Ivy_ManConst1(pNew);
    Ivy_ManForEachPi( p, pObj, i )
        pObj->pEquiv = Ivy_ObjCreatePi(pNew);
    // create the fake PIs for latches
    Ivy_ManForEachNodeVec( p, vLatches, pObj, i )
        pObj->pEquiv = Ivy_ObjCreatePi(pNew);
    // duplicate internal nodes
    Ivy_ManForEachNodeVec( p, vNodes, pObj, i )
        pObj->pEquiv = Ivy_And( pNew, Ivy_ObjChild0Equiv(pObj), Ivy_ObjChild1Equiv(pObj) );
    // add the POs
    Ivy_ManForEachPo( p, pObj, i )
        Ivy_ObjCreatePo( pNew, Ivy_ObjChild0Equiv(pObj) );
    // transform additional PI nodes into latches and connect them
    Ivy_ManForEachNodeVec( p, vLatches, pObj, i )
    {
        assert( !Ivy_ObjFaninC0(pObj) );
        pObj->pEquiv->Type = IVY_LATCH;
        pObj->pEquiv->Init = pObj->Init;
        Ivy_ObjConnect( pNew, pObj->pEquiv, Ivy_ObjChild0Equiv(pObj), NULL );
    }
    // shrink the arrays
    Vec_PtrShrink( pNew->vPis, Ivy_ManPiNum(p) );
    // update the counters of different objects
    pNew->nObjs[IVY_PI] -= Ivy_ManLatchNum(p);
    pNew->nObjs[IVY_LATCH] += Ivy_ManLatchNum(p);
    // free arrays
    Vec_IntFree( vNodes );
    Vec_IntFree( vLatches );
    // make sure structural hashing did not change anything
    assert( Ivy_ManNodeNum(p)  == Ivy_ManNodeNum(pNew) );
    assert( Ivy_ManLatchNum(p) == Ivy_ManLatchNum(pNew) );
    // check the resulting network
    if ( !Ivy_ManCheck(pNew) )
        printf( "Ivy_ManMakeSeq(): The check has failed.\n" );
    return pNew;
}

/**Function*************************************************************

  Synopsis    [Stops the AIG manager.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Ivy_ManStop( Ivy_Man_t * p )
{
    if ( p->time1 ) { PRT( "Update lev  ", p->time1 ); }
    if ( p->time2 ) { PRT( "Update levR ", p->time2 ); }
//    Ivy_TableProfile( p );
//    if ( p->vFanouts )  Ivy_ManStopFanout( p );
    if ( p->vChunks )   Ivy_ManStopMemory( p );
    if ( p->vRequired ) Vec_IntFree( p->vRequired );
    if ( p->vPis )      Vec_PtrFree( p->vPis );
    if ( p->vPos )      Vec_PtrFree( p->vPos );
    if ( p->vBufs )     Vec_PtrFree( p->vBufs );
    if ( p->vObjs )     Vec_PtrFree( p->vObjs );
    free( p->pTable );
    free( p );
}

/**Function*************************************************************

  Synopsis    [Returns the number of dangling nodes removed.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Ivy_ManCleanup( Ivy_Man_t * p )
{
    Ivy_Obj_t * pNode;
    int i, nNodesOld;
    nNodesOld = Ivy_ManNodeNum(p);
    Ivy_ManForEachObj( p, pNode, i )
        if ( Ivy_ObjIsNode(pNode) || Ivy_ObjIsLatch(pNode) || Ivy_ObjIsBuf(pNode) )
            if ( Ivy_ObjRefs(pNode) == 0 )
                Ivy_ObjDelete_rec( p, pNode, 1 );
    return nNodesOld - Ivy_ManNodeNum(p);
}

/**Function*************************************************************

  Synopsis    [Returns the number of dangling nodes removed.]

  Description []
               
  SideEffects []

  SeeAlso     [] 

***********************************************************************/
int Ivy_ManPropagateBuffers( Ivy_Man_t * p, int fUpdateLevel )
{
    Ivy_Obj_t * pNode;
    int nSteps;
    for ( nSteps = 0; Vec_PtrSize(p->vBufs) > 0; nSteps++ )
    {
        pNode = Vec_PtrEntryLast(p->vBufs);
        while ( Ivy_ObjIsBuf(pNode) )
            pNode = Ivy_ObjReadFirstFanout( p, pNode );
        Ivy_NodeFixBufferFanins( p, pNode, fUpdateLevel );
    }
//    printf( "Number of steps = %d\n", nSteps );
    return nSteps;
}

/**Function*************************************************************

  Synopsis    [Stops the AIG manager.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Ivy_ManPrintStats( Ivy_Man_t * p )
{
    printf( "PI/PO = %d/%d ", Ivy_ManPiNum(p), Ivy_ManPoNum(p) );
    printf( "A = %7d. ",       Ivy_ManAndNum(p) );
    printf( "L = %5d. ",       Ivy_ManLatchNum(p) );
//    printf( "X = %d. ",       Ivy_ManExorNum(p) );
//    printf( "B = %3d. ",       Ivy_ManBufNum(p) );
    printf( "MaxID = %7d. ",   Ivy_ManObjIdMax(p) );
//    printf( "Cre = %d. ",     p->nCreated );
//    printf( "Del = %d. ",     p->nDeleted );
    printf( "Lev = %3d. ",     Ivy_ManLatchNum(p)? -1 : Ivy_ManLevels(p) );
    printf( "\n" );
}

/**Function*************************************************************

  Synopsis    [Converts a combinational AIG manager into a sequential one.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Ivy_ManMakeSeq( Ivy_Man_t * p, int nLatches, int * pInits )
{
    Ivy_Obj_t * pObj, * pLatch;
    Ivy_Init_t Init;
    int i;
    if ( nLatches == 0 )
        return;
    assert( nLatches < Ivy_ManPiNum(p) && nLatches < Ivy_ManPoNum(p) );
    assert( Ivy_ManPiNum(p) == Vec_PtrSize(p->vPis) );
    assert( Ivy_ManPoNum(p) == Vec_PtrSize(p->vPos) );
    assert( Vec_PtrSize( p->vBufs ) == 0 );
    // create fanouts
    if ( p->fFanout == 0 )
        Ivy_ManStartFanout( p );
    // collect the POs to be converted into latches
    for ( i = 0; i < nLatches; i++ )
    {
        // get the latch value
        Init = pInits? pInits[i] : IVY_INIT_0;
        // create latch
        pObj = Ivy_ManPo( p, Ivy_ManPoNum(p) - nLatches + i );
        pLatch = Ivy_Latch( p, Ivy_ObjChild0(pObj), Init );
        Ivy_ObjDisconnect( p, pObj );
        // recycle the old PO object
        Vec_PtrWriteEntry( p->vObjs, pObj->Id, NULL );
        Ivy_ManRecycleMemory( p, pObj );
        // convert the corresponding PI to a buffer and connect it to the latch
        pObj = Ivy_ManPi( p, Ivy_ManPiNum(p) - nLatches + i );
        pObj->Type = IVY_BUF;
        Ivy_ObjConnect( p, pObj, pLatch, NULL );
        // save the buffer
        Vec_PtrPush( p->vBufs, pObj );
    }
    // shrink the arrays
    Vec_PtrShrink( p->vPis, Ivy_ManPiNum(p) - nLatches );
    Vec_PtrShrink( p->vPos, Ivy_ManPoNum(p) - nLatches );
    // update the counters of different objects
    p->nObjs[IVY_PI] -= nLatches;
    p->nObjs[IVY_PO] -= nLatches;
    p->nObjs[IVY_BUF] += nLatches;
    p->nDeleted -= 2 * nLatches;
    // remove dangling nodes
    Ivy_ManCleanup(p);
/* 
    // check for dangling nodes
    Ivy_ManForEachObj( p, pObj, i )
        if ( !Ivy_ObjIsPi(pObj) && !Ivy_ObjIsPo(pObj) && !Ivy_ObjIsConst1(pObj) )
        {
            assert( Ivy_ObjRefs(pObj) > 0 );
            assert( Ivy_ObjRefs(pObj) == Ivy_ObjFanoutNum(p, pObj) );
        }
*/
    // perform hashing by propagating the buffers
    Ivy_ManPropagateBuffers( p, 0 );
    // fix the levels
    Ivy_ManResetLevels( p );
    // check the resulting network
    if ( !Ivy_ManCheck(p) )
        printf( "Ivy_ManMakeSeq(): The check has failed.\n" );
}

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


