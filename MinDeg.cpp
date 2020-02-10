//-----------------------------------------------
// Class for graph manipulations and tools for sparse matrix 
// Reordering. The Reordering is performed to minimize the
// number of nonzeros in the factorized matrix.
// Algorithms are based on the reordering algorithms 
// proposed by Alan George and  JWH Liu.
// See https://epubs.siam.org/doi/abs/10.1137/1031001
//
// By: Mehdi Paak
// TODO: major refactoring needed
//----------------------------------------------

#include "stdafx.h"
#pragma once

#include<iostream>
#include<cassert>
#include<sys/timeb.h>
#include "MinDeg.h"

//=======================================================
// Graph
// 
//=======================================================
// Constructor
 REORD_QGraph:: REORD_QGraph(
  const int i_nNumV, 
  const int i_nNumE,
  Indexing i_Style)
{

  m_nNumV = i_nNumV;
  m_nNumE = i_nNumE;
  m_nStyle = i_Style;

  SetDataSpace(m_nNumV, m_nNumE);
}

//========================================
// Copy constructor

 REORD_QGraph:: REORD_QGraph(const  REORD_QGraph &i_sGraph)
{
  // Check input.
  assert(i_sGraph.IsEmpty() == false);

  m_nNumV = i_sGraph.GetNumV();
  m_nNumE = i_sGraph.GetNumE();
  m_nStyle = i_sGraph.GetStyle();

  SetDataSpace(m_nNumV, m_nNumE);
  
  ::memcpy(m_panXAdj, i_sGraph.GetXAdj(), sizeof(*m_panXAdj) * (m_nNumV + 1));
  ::memcpy(m_panAdj, i_sGraph.GetAdj(), sizeof(*m_panAdj) * m_nNumE);
}

//========================================
// Destructor
 REORD_QGraph::~ REORD_QGraph()
{
  printf("Deleting QGraph\n");
  
  if(m_nStyle == F_Style)
  {
    m_panXAdj++;
    m_panAdj++;
  }

  if(m_panXAdj)
  {
    delete [] m_panXAdj;
    m_panXAdj = 0;
  }

  if(m_panAdj)
  {
    delete [] m_panAdj;
    m_panAdj = 0;
  }

}

//========================================
// Set Data
void  REORD_QGraph::SetDataSpace(
  const int i_nNumV, 
  const int i_nNumE)
{
  // Check input.
  assert(i_nNumV > 0 && i_nNumE >= 0);

  try
  {  
    m_panXAdj = new int [i_nNumV + 1];
    
    if(i_nNumE > 0)
      m_panAdj = new int [i_nNumE];
  }

  catch(...)
  {
    assert(false);
    delete [] m_panXAdj;
    
    if(i_nNumE > 0)
      delete [] m_panAdj;
    throw;
  }
}

//=====================================
// Print Graph

void  REORD_QGraph::PrintGraph() const
{
  
  assert(!IsEmpty());

  // Set the start and end index appropriate to C or F style. 
  int nStartInd, nStopInd;
  if(m_nStyle == C_Style)
  {
    nStartInd = 0; nStopInd = m_nNumV;
  }
  else if(m_nStyle == F_Style)
  {
    nStartInd = 1; nStopInd = m_nNumV + 1;
  }

  printf("\nNode -- Graph\n");
  for(int i = nStartInd; i < nStopInd; i++)
  {
    printf("%-3d --> ", i);
    for(int j = m_panXAdj[i]; j < m_panXAdj[i + 1]; j++)
      printf("%-2d ", m_panAdj[j]);
    printf("\n");
  }
 
/*  printf("\nNode-- Deg -- Marker -- QSize -- QLink--\n");
  for(int i = nStartInd; i < nStopInd; i++)
  {
    printf("%-2d     %-2d       %-2d        %-2d       %-2d", i, m_panDeg[i], m_panMarker[i], m_panQSize[i], m_panQLink[i]);
    printf("\n");
  }*/
  printf("\n======================\n");
}


void  REORD_QGraph::PrintToFile(const char *i_pcFileName)
{
  // Check input.
  assert(i_pcFileName);

  // Set the start and end index appropriate to C or F style. 
  const int nNNZ = this->GetNumE();
  int nStart, nStopPtr, nStopInd;
  if(m_nStyle == C_Style)
  {
    nStart = 0; 
    nStopPtr = m_nNumV + 1;
    nStopInd = nNNZ;
  }
  else if(m_nStyle == F_Style)
  {
    nStart = 1; 
    nStopPtr = m_nNumV + 2;
    nStopInd = nNNZ + 1;
  }

  FILE *pFile;
  if(fopen_s(&pFile, i_pcFileName,"w"))
    fprintf(stderr, "Error opening File %s \n", i_pcFileName);

  if (m_panXAdj)
  {
   for (int i = nStart; i < nStopPtr; i++)
     fprintf(pFile,"%d ",m_panXAdj[i]); 
   fprintf(pFile,"\n");
  }

  if (m_panAdj)
  {
   for (int i = nStart; i <  nStopInd; i++)
     fprintf(pFile,"%d ", m_panAdj[i]); 
   fprintf(pFile,"\n");
  }

/*  if (m_padVal)
  {
   for (int i = nStart; i <  nStopInd; i++)
     _ftprintf_s(pFile,_T("%-22.16f "),m_padValue[i]); 
  } */

  fclose(pFile);
}
// Convert to F_Style indexing
void  REORD_QGraph::ConvertToFort()
{
  // Check input.
  assert(!IsEmpty() && GetStyle() == C_Style);

  // Add 1 to XAdj, Adj: i.e., vertex 0 becomes 1 and so on. 0 -> 1.

  for(int i = 0; i < m_nNumV + 1; i++)
  {
    m_panXAdj[i]++;
  }

  for(int i = 0; i < m_nNumE; i++)
  {
    m_panAdj[i]++;
  }

  // Shift the pointer one back so that Arr[1] (F style) points to Arr[0] (C style).
  m_panXAdj--;
  m_panAdj--;

  // Set Style
  SetStyleTo(F_Style);
}

// Convert to C_Style indexing
void  REORD_QGraph::ConvertToC()
{
  // Check input.
  assert(!IsEmpty() && GetStyle() == F_Style);

  // Subtract 1 from XAdj, Adj: i.e., vertex 0 becomes 1 and so on. 0 -> 1.

  for(int i = 1; i <= m_nNumV + 1; i++)
  {
    m_panXAdj[i]--;
  }

  for(int i = 1; i <= m_nNumE; i++)
  {
    m_panAdj[i]--;
  }

  // Shift the pointer one back so that Arr[1] (F style) points to Arr[0] (C style).
  m_panXAdj++;
  m_panAdj++;

  // Set Style
  SetStyleTo(C_Style);
}

// Read graph from file.
void  REORD_QGraph::ReadGrphFromFile(
   REORD_QGraph *i_psQGraph,
  const char *i_pacFileName)
{
  //Check input.
  assert(i_psQGraph && !i_psQGraph->IsEmpty() && i_pacFileName);

  const int nNumV = i_psQGraph->GetNumV();
  const int nNumE = i_psQGraph->GetNumE();

  int *__restrict panXAdj = i_psQGraph->GetXAdj();
  int *__restrict panAdj = i_psQGraph->GetAdj(); 

  FILE * pFile;
  int nError = fopen_s(&pFile,i_pacFileName,"r");
  if(nError != 0)
  {
    printf("\nError opening File\n");
    return;
  }

  for(int i = 0; i < (nNumV + 1); i++)
    fscanf_s(pFile, "%d", &panXAdj[i]);

  for(int i = 0; i < nNumE; i++)
    fscanf_s(pFile, "%d", &panAdj[i]);

  fclose(pFile);
}

// Convert unordered graph with diagonals to QGraph format rewrite into file.
// diagonals should be removed to avoid self-edge.
void  REORD_QGraph::ChangeToQGrph(
  const int &i_nNumV,
  const int &i_nNumE,
  const char *i_pacInFileName,
  const char *i_pacOutFileName)
{
  //Check input.
  assert(i_nNumV > 0 && i_pacInFileName && i_pacOutFileName);

  FILE * pInFile;
  int nError = fopen_s(&pInFile,i_pacInFileName,"r");
  if(nError != 0)
  {
    printf("\nError opening Input File\n");
    return;
  }

  FILE * pOutFile;
  nError = fopen_s(&pOutFile,i_pacOutFileName,"w");
  if(nError != 0)
  {
    printf("\nError opening Output File\n");
    return;
  }

  int *__restrict panPtr = new(std::nothrow) int [i_nNumV + 1];
  if(!panPtr)
  {
  assert(false);
  return;
  }
  
  for(int i = 0; i < (i_nNumV + 1); i++)
  {
    int nPtrVal;
    fscanf_s(pInFile, "%d", &nPtrVal);
    panPtr[i] = nPtrVal;
    fprintf_s(pOutFile,"%d ",nPtrVal - i);
    
  }
  
  fprintf_s(pOutFile,"\n");

  for(int i = 0; i < i_nNumV + 1; i++)
  {
    for(int j = panPtr[i]; j < panPtr[i + 1]; j++)
    {
      int nIndVal;
      fscanf_s(pInFile, "%d ", &nIndVal);
      if(nIndVal == i)
        continue;
      fprintf_s(pOutFile,"%d ", nIndVal);
    }
  }
  fclose(pInFile);
  fclose(pOutFile);
  delete [] panPtr;
}

//++++++++++++++++++++++++++
// Sparse Mat
//
//
//++++++++++++++++++++++++++

// Constructor.
SparseMat::SparseMat(
  const int i_nNumV, 
  const int i_nNumE, 
  Indexing i_Style):  REORD_QGraph(i_nNumV,i_nNumE,i_Style)
{
  SetDataSpaceSparseMat(i_nNumE);
}

// Copy constructor.
SparseMat::SparseMat(const SparseMat &i_psSparseMat)
{
  // Check input.
  assert(i_psSparseMat.IsEmpty() == false);

  m_nNumV = i_psSparseMat.GetNumV();
  m_nNumE = i_psSparseMat.GetNumE();
  m_nStyle = i_psSparseMat.GetStyle();

  SetDataSpace(m_nNumV, m_nNumE);
  SetDataSpaceSparseMat(m_nNumE);

  ::memcpy(m_panXAdj, i_psSparseMat.GetXAdj(), sizeof(*m_panXAdj) * (m_nNumV + 1));
  ::memcpy(m_panAdj,  i_psSparseMat.GetAdj(), sizeof(*m_panAdj) * m_nNumE);
  ::memcpy(m_padVal,  i_psSparseMat.GetVal(), sizeof(*m_padVal) * m_nNumE);

}

// Destructor.
SparseMat::~SparseMat()
{

  printf("Deleting SparseMat\n");
  
  if(m_nStyle == F_Style)
  {
    m_padVal++;
  }

  if(m_padVal)
  {
    delete [] m_padVal;
    m_padVal = 0;
  }

  // Base class destructor automatically called.

}

// Set Data.
void SparseMat::SetDataSpaceSparseMat(const int i_nNumE)
{

  // Check input.
  assert(i_nNumE >= 0);

  try
  {  
    m_padVal = new double [i_nNumE];
  }

  catch(...)
  {
    assert(false);
    delete [] m_padVal;
    throw;
  }
}

//================================================================
// Form graph adjaceny from matrix upper triangular part ordered.
// 
//================================================================
void SparseMat::FormGraphAdj(
  const SparseMat *i_psU,
   REORD_QGraph *o_psG,
  int &o_nMaxRowNNZ)
{
  // Check input.
  assert(i_psU && o_psG && (i_psU->GetStyle() == C_Style) && (o_psG->GetStyle() == C_Style));

  const int nNumV = i_psU->GetNumV();
  const int *panPtrU = i_psU->GetXAdj();
  const int *panIndxU = i_psU->GetAdj();

  const int nNumGE = o_psG->GetNumE();
  int *panPtrG = o_psG->GetXAdj();
  int *panIndxG = o_psG->GetAdj();
  

  // Initiallization.
  ::memset(panPtrG, 0, sizeof(*panPtrG) * (nNumV+1));
  ::memset(panIndxG, 0, sizeof(*panIndxG) * nNumGE);

  // Adjacency list row sum in PtrG.
  o_nMaxRowNNZ = 0;
  for(int i = 0; i < nNumV; i++)
  {
    const int nPtrUStrt = panPtrU[i] + 1; // Remove diagonal from count (no self-edge).
    const int nPtrUStop = panPtrU[i + 1];
    const int nSum1 = nPtrUStop - nPtrUStrt;

    for(int j = nPtrUStrt; j < nPtrUStop; j++)
    {
      panPtrG[panIndxU[j]]++;
    }
       
    panPtrG[i] += nSum1;
    if(panPtrG[i] > o_nMaxRowNNZ)
      o_nMaxRowNNZ = panPtrG[i];
  }

  // Cumulative sum PtrG.
  int nS0 = panPtrG[0];
  int nS1 = nS0;
  panPtrG[0] = 0;
  for(int i = 1; i < nNumV + 1; i++)
  { 
    nS0 = panPtrG[i];
    panPtrG[i] = panPtrG[i - 1] + nS1; 
    nS1 = nS0;
  }

  // Go through Indices of U and place appropriately in IndxG; exclude diagonal.
  
  int *panWrk = new int[nNumV];
  ::memcpy(panWrk,panPtrG,sizeof(*panWrk) * nNumV);

  for(int i = 0; i < nNumV; i++)
  {
    for(int j = panPtrU[i] + 1; j < panPtrU[i + 1]; j++)
    {
      const int nIndxU = panIndxU[j];
      panIndxG[panWrk[i]++] = nIndxU; // copy Row i form U to G.
      panIndxG[panWrk[nIndxU]++] = i; // Update other Rows in G  which are symmtric in Row i of U.
    }
       
  }

  // Free memory.
  delete [] panWrk;
}

//================================================================
// Get Diagonal from the sorted upper triangular matrix.
// 
//================================================================
void SparseMat::GetDiagonal(
  const SparseMat *i_psU, 
  const MatForm i_nMatForm, 
  double *o_padDiag)
{
  // Check input.
  assert(i_psU && i_nMatForm == nUpper && o_padDiag);

  const int nNumEq = i_psU->GetNumV();
  const int *panPtr = i_psU->GetXAdj();
  const double *padVal = i_psU->GetVal();

  for(int i = 0; i < nNumEq; i++)
  {
    const int j = panPtr[i];
    o_padDiag[i] = padVal[j];
  }
}

//================================================================
// Insert A into L (NNZL>>NNZA) based on permutation vec, 
// Diagonals are stored separately. A = LL^T
// 
//================================================================
void SparseMat::InsertAToL(
  const SparseMat *i_psU, 
  const MatForm i_nMatForm, 
  const int *i_panPerm, 
  const int *i_panInvPerm, 
  double *o_padLnz,
  double *o_padDiag)
{
  // Check input.
  assert(i_psU && i_nMatForm == nUpper && o_padDiag);





}


//==========================
//
// MinDegData
//==========================

// Constructor.
  REORD_MinDegData_s:: REORD_MinDegData_s(const int i_nNumV)
{
  // Check input.
  assert(i_nNumV > 0);

  try
  {
    m_panDeg = new int [i_nNumV];
    m_panMarker = new int [i_nNumV];
    m_panQLink = new int [i_nNumV];
    m_panQSize = new int [i_nNumV];

    m_panRchSet = new int [i_nNumV];
    m_panNbrhdSet = new int [i_nNumV];
    m_nRchSize = 0;
    m_nNbrhdSize = 0;
    m_nRchTotNodes = 0;
    m_nNumElimSupNodes = 0;

    m_nNumV = i_nNumV;
    m_nStyle = C_Style;
  }
  catch(...)
  {
    delete [] m_panDeg;
    delete [] m_panMarker;
    delete [] m_panQLink;
    delete [] m_panQSize;
    
    delete [] m_panRchSet;
    delete [] m_panNbrhdSet;
    assert(false);
  }

}

// Destructor.
 REORD_MinDegData_s::~ REORD_MinDegData_s()
{ 
  printf("\nDeleting MinDegData\n");

  // Adjusting pointer based on indexing style.
  if(m_nStyle == F_Style)
  {
    m_panDeg++;
    m_panMarker++;
    m_panQLink++;
    m_panQSize++;

    m_panRchSet++;
    m_panNbrhdSet++;
  }

  if(m_panDeg)
  {
    delete [] m_panDeg;
    m_panDeg = 0;
  }

  if(m_panMarker)
  {
    delete [] m_panMarker;
    m_panMarker = 0;
  }

  if(m_panQLink)
  {
    delete [] m_panQLink;
    m_panQLink = 0;
  }

  if(m_panQSize)
  {
    delete [] m_panQSize;
    m_panQSize = 0;
  }

  if(m_panRchSet)
  {
    delete [] m_panRchSet;
    m_panRchSet = 0;
  }

  if(m_panNbrhdSet)
  {
    delete [] m_panNbrhdSet;
    m_panNbrhdSet = 0;
  }
}

void  REORD_MinDegData_s::ConvertToFort()
{
  // Check input.
  assert(!IsEmpty() && GetStyle() == C_Style);
  
  // Shift Pointer.
  m_panDeg--;
  m_panMarker--;
  m_panQLink--;
  m_panQSize--;

  m_panRchSet--;
  m_panNbrhdSet--;
  SetStyleTo(F_Style);
}

void  REORD_MinDegData_s::ConvertToC()
{
  // Check input.
  assert(!IsEmpty() && GetStyle() == F_Style);
  
  // Shift Pointer.
  m_panDeg++;
  m_panMarker++;
  m_panQLink++;
  m_panQSize++;

  m_panRchSet++;
  m_panNbrhdSet++;
  SetStyleTo(C_Style);
}

// Print Vec
void  REORD_MinDegData_s::PrintVec(
  const int &i_nLen, 
  const int *i_panVec, 
  const Indexing i_Style)
{
  // Check data.
  assert(i_panVec);
   if(i_nLen <= 0) return;

  // Set the start and end index appropriate to C or F style. 
  int nStartInd, nStopInd;
  if(i_Style == C_Style)
  {
    nStartInd = 0; nStopInd = i_nLen;
  }
  else if(i_Style == F_Style)
  {
    nStartInd = 1; nStopInd = i_nLen + 1;
  }

  printf("\n------------\n");

  for(int i = nStartInd; i < nStopInd; i++)
    printf("%3d -- %4d\n", i, i_panVec[i]);

}

//=========================
// Initialize data for MD
//
//=========================
void  REORD_SDSMinDeg::InitDataForMD(
   REORD_QGraph *io_psQGraph, 
   REORD_MinDegData_s *o_psMinDegData,
  int * __restrict o_panPerm,
  int * __restrict o_panInvPerm,
  int &o_nFirstMinDeg)
{
  // Check input.
  assert(io_psQGraph && o_psMinDegData && io_psQGraph->GetStyle() == F_Style && o_psMinDegData->GetStyle() == F_Style);

  const int nNumV = io_psQGraph->GetNumV(); // Number of vertices/nodes.
  

  // Quotient Graph.
  const int * __restrict panXAdj = io_psQGraph->GetXAdj();
  assert(panXAdj);
  int * __restrict panMarker = o_psMinDegData->GetMarker();
  assert(panMarker);
  int * __restrict panQLink = o_psMinDegData->GetQLink();
  assert(panQLink);
  int * __restrict panQSize = o_psMinDegData->GetQSize();
  assert(panQSize);
  int * __restrict panDeg = o_psMinDegData->GetDeg();
  assert(panDeg);

  // MinDegData
  int * __restrict panRchSet = o_psMinDegData->GetRchSet();
  assert(panRchSet);
  int * __restrict panNbrhdSet = o_psMinDegData->GetNbrhdSet();
  assert(panNbrhdSet);

  
  o_nFirstMinDeg = nNumV;

  // Initialization loop, All Parameters must have F_Style.
  for(int i = 1; i <= nNumV; i++)
  {
    o_panPerm[i] = i;
    o_panInvPerm[i] = i;
    panMarker[i] = 0;
    panQSize[i] = 1;
    panQLink[i] = 0;
    
    panRchSet[i] = 0;
    panNbrhdSet[i] = 0;

    const int nDeg = panXAdj[i + 1] - panXAdj[i];
    panDeg[i] = nDeg;
    if(o_nFirstMinDeg > nDeg)
      o_nFirstMinDeg = nDeg;
  }

}

//===========================
// Find Min Deg permutation.
//===========================
void  REORD_SDSMinDeg::FindMinDegPerm(
   REORD_QGraph *io_psQGraph,
  int * __restrict o_panPerm,
  int * __restrict o_panInvPerm,
  SDS_INT & o_nNumSubs)
{
  // Check input.
  assert(io_psQGraph && !io_psQGraph->IsEmpty() && o_panPerm && o_panInvPerm);

  const int nNumV = io_psQGraph->GetNumV();

  // Define MinDegData.
   REORD_MinDegData_s sMinDegData(nNumV); // TODO: check Free Mem.

  // Convert to F_Style.
  io_psQGraph->ConvertToFort();
  sMinDegData.ConvertToFort();
  
  o_panPerm--;
  o_panInvPerm--;

  // Initialize.
  int nMinDeg = 0;   // Used for threshold search George & Liu.
  int nThrshld = 0;  // Used for threshold search George & Liu.
  o_nNumSubs = 0;    // Number of compressed subscripts to be used for compressed index format.

   REORD_SDSMinDeg::InitDataForMD(
    io_psQGraph,
    &sMinDegData,
    o_panPerm,
    o_panInvPerm,
    nMinDeg);

  //
  int * __restrict panMarker = sMinDegData.GetMarker();
  assert(panMarker);
  int * __restrict panDeg = sMinDegData.GetDeg();
  assert(panDeg);
  int * __restrict panQLink = sMinDegData.GetQLink();
  assert(panQLink);
  int * __restrict panQSize = sMinDegData.GetQSize();

  // Main loop: find node of min deg, eliminate, update deg, update quotient graph, repeat.
  int nCurrNode = 0;     // Current node.
  int nCurrNum = 0;      // Current number.
  int nSrchPos = 0;      // Search start position.
  int nNNZ = 0;          // Number of nonzeros in U factor. 

  bool bSrchFlag1 = true; 
  bool bSrchFlag2 = true;
  
  while(nCurrNum < nNumV)
  {
    
    // Find Node/XNode with min deg.
    while(bSrchFlag1)
    {
      
      if(bSrchFlag2)
      {
        nSrchPos = 1;
        nThrshld = nMinDeg;
        nMinDeg = nNumV;
      }

      const int nNextPosinPerm = nCurrNum + 1;
    
      // In case XNodes eliminated the number of nodes in XNode is taken into account.
      if(nNextPosinPerm > nSrchPos) 
        nSrchPos = nNextPosinPerm;

      for(int j = nSrchPos; j <= nNumV; j++)
      {
        nCurrNode = o_panPerm[j];
        if (panMarker[nCurrNode] < 0) 
          continue;
        const int nDeg = panDeg[nCurrNode];
        if(nDeg <= nThrshld)
        {
          bSrchFlag1 = false;  // Min Deg node found, break out of while loop by Flag 1 = false
          bSrchFlag2 = false;  // Next round do not reset search vars by Flag2 = false.
          nSrchPos = j;        // Next Srch atart position if no XNode or modified Srch param by nodes in RchSet.
          o_nNumSubs += nDeg;     // TODO: export info.

          const int nQsz = panQSize[nCurrNode];
          nNNZ += nDeg * nQsz - nQsz * (nQsz - 1)/2; // TODO: export info.
          
          break;
        }
        else if(nMinDeg > nDeg)
        {
          nMinDeg = nDeg;
        }
      }

      if(bSrchFlag1) 
        bSrchFlag2 = true; // In case no Min Deg node found within Srch--NumV, Flag2=true enables resetting to Srch = 1 etc and searching from the beginning.

    } // End Srch While.

    bSrchFlag1 = true;  // Enable the search loop for the next round.

   // io_psQGraph->PrintGraph();

    // Find ReachSet(Node, S).
    panMarker[nCurrNode] = 1;  // Set to 1 to avoid including itself in RchSet, will be reset back to 0 before QGUpdate.
    FindReachSet(
      nCurrNode,
      io_psQGraph,
      &sMinDegData);

    const int nRchSize = sMinDegData.GetRchSize();
    const int nNbrhdSize = sMinDegData.GetNbrhdSize();
    const int * __restrict panRchSet = sMinDegData.GetRchSet();

    // Eliminate CurrNode/CurrXNode, go through Qlink for XNodes.
    int nXnode = nCurrNode;
    while(nXnode > 0)
    {
      nCurrNum++;

      const int nPosInPerm = o_panInvPerm[nXnode];
      const int nPosInInv = o_panPerm[nCurrNum]; 

      o_panPerm[nPosInPerm] = nPosInInv;
      o_panInvPerm[nPosInInv] = nPosInPerm;

      o_panPerm[nCurrNum] = nXnode;
      o_panInvPerm[nXnode] = nCurrNum;

      panDeg[nXnode] = -1; // Eliminated.
      nXnode = panQLink[nXnode];
    }


    if(nRchSize > 0)
    {

      // Update degree of Nodes in RchSet and find possible XNodes.
      UpdateDeg(
        io_psQGraph,
        &sMinDegData);

      // Reset RchSet markers from 2 (in update) to 0 except for non-representative nodes in an XNode with marker = -1. 
      // Check if any of RchSet node has min deg, if so set threshold, SrchPos etc.

      panMarker[nCurrNode] = 0; // Current eliminated supernode representative.
      for(int i = 1; i <= nRchSize; i++)
      {
        const int nRchSetNode = panRchSet[i];
        if(panMarker[nRchSetNode] < 0)
          continue;
        panMarker[nRchSetNode] = 0;
        
        const int nDeg = panDeg[nRchSetNode];
        if(nDeg < nMinDeg) 
          nMinDeg = nDeg;
        if(nDeg <= nThrshld)
        {
          nMinDeg = nThrshld;
          nThrshld = nDeg;
          nSrchPos = o_panInvPerm[nRchSetNode];
        }
      }

      // Update Quotient Graph.
      if(nNbrhdSize > 0)
        TransformQuotientGraph(
        nCurrNode,
        io_psQGraph,
        &sMinDegData);
    }

  } // End While over nodes.

  // Convert Perm and InvPerm back to C_Style indexing, Indices are shifted back to start from 0.
  // Shift Values.
  for(int i = 1; i <= nNumV; i++)
  {
    o_panPerm[i]--;
    o_panInvPerm[i]--;
  }
  // Shift Pointers.
  o_panPerm++; 
  o_panInvPerm++;

}

//===========================
// Find Min Deg permutation.
//===========================
void  REORD_SDSMinDeg::FindMinDegPerm_2(
   REORD_QGraph *io_psQGraph,
  int * __restrict o_panPerm)
{
  // Check input.
  assert(io_psQGraph && !io_psQGraph->IsEmpty() && o_panPerm);

  const int nNumV = io_psQGraph->GetNumV();

  // Define MinDegData.
   REORD_MinDegData_s sMinDegData(nNumV); // TODO: check Free Mem.
  int * __restrict panInvPerm = new(std::nothrow) int [nNumV];
  assert(panInvPerm);

  // Convert to F_Style.
  io_psQGraph->ConvertToFort();
  sMinDegData.ConvertToFort();
  
  o_panPerm--;
  panInvPerm--;

  // Initialize.
  int nMinDeg = 0;   // Used for threshold search George & Liu.
  int nThrshld = 0;  // Used for threshold search George & Liu.

   REORD_SDSMinDeg::InitDataForMD(
    io_psQGraph,
    &sMinDegData,
    o_panPerm,
    panInvPerm,
    nMinDeg);

  //
  int * __restrict panMarker = sMinDegData.GetMarker();
  assert(panMarker);
  int * __restrict panDeg = sMinDegData.GetDeg();
  assert(panDeg);
  int * __restrict panQLink = sMinDegData.GetQLink();
  assert(panQLink);
  int * __restrict panQSize = sMinDegData.GetQSize();

  // Main loop: find node of min deg, eliminate, update deg, update quotient graph, repeat.
  int nCurrNode = 0;     // Current node.
  int nCurrNum = 0;      // Current number.
  int nSrchPos = 0;      // Search start position.
  int nNumSub = 0;       // Number of substitutions.
  int nNNZ = 0;          // Number of nonzeros in U factor. 

  bool bSrchFlag1 = true; 
  bool bSrchFlag2 = true;
  
  while(nCurrNum < nNumV)
  {
    
    // Find Node/XNode with min deg.
    while(bSrchFlag1)
    {
      
      if(bSrchFlag2)
      {
        nSrchPos = 1;
        nThrshld = nMinDeg;
        nMinDeg = nNumV;
      }

      const int nNextPosinPerm = nCurrNum + 1;
    
      // In case XNodes eliminated the number of nodes in XNode is taken into account.
      if(nNextPosinPerm > nSrchPos) 
        nSrchPos = nNextPosinPerm;

      for(int j = nSrchPos; j <= nNumV; j++)
      {
        nCurrNode = o_panPerm[j];
        if (panMarker[nCurrNode] < 0) 
          continue;
        const int nDeg = panDeg[nCurrNode];
        if(nDeg <= nThrshld)
        {
          bSrchFlag1 = false;  // Min Deg node found, break out of while loop by Flag 1 = false
          bSrchFlag2 = false;  // Next round do not reset search vars by Flag2 = false.
          nSrchPos = j;        // Next Srch atart position if no XNode or modified Srch param by nodes in RchSet.
          nNumSub += nDeg;     // TODO: export info.

          const int nQsz = panQSize[nCurrNode];
          nNNZ += nDeg * nQsz - nQsz * (nQsz - 1)/2; // TODO: export info.
          
          break;
        }
        else if(nMinDeg > nDeg)
        {
          nMinDeg = nDeg;
        }
      }

      if(bSrchFlag1) 
        bSrchFlag2 = true; // In case no Min Deg node found within Srch--NumV, Flag2=true enables resetting to Srch = 1 etc and searching from the beginning.

    } // End Srch While.

    bSrchFlag1 = true;  // Enable the search loop for the next round.

   // io_psQGraph->PrintGraph();

    // Find ReachSet(Node, S).
    panMarker[nCurrNode] = 1;  // Set to 1 to avoid including itself in RchSet, will be reset back to 0 before QGUpdate.
    FindReachSet(
      nCurrNode,
      io_psQGraph,
      &sMinDegData);

    const int nRchSize = sMinDegData.GetRchSize();
    const int nNbrhdSize = sMinDegData.GetNbrhdSize();
    const int * __restrict panRchSet = sMinDegData.GetRchSet();

    // Eliminate CurrNode/CurrXNode, go through Qlink for XNodes.
    int nXnode = nCurrNode;
    while(nXnode > 0)
    {
      nCurrNum++;

      const int nPosInPerm = panInvPerm[nXnode];
      const int nPosInInv = o_panPerm[nCurrNum]; 

      o_panPerm[nPosInPerm] = nPosInInv;
      panInvPerm[nPosInInv] = nPosInPerm;

      o_panPerm[nCurrNum] = nXnode;
      panInvPerm[nXnode] = nCurrNum;

      panDeg[nXnode] = -1; // Eliminated.
      nXnode = panQLink[nXnode];
    }


    if(nRchSize > 0)
    {

      // Update degree of Nodes in RchSet and find possible XNodes.
      UpdateDeg(
        io_psQGraph,
        &sMinDegData);

      // Reset RchSet markers from 2 (in update) to 0 except for non-representative nodes in an XNode with marker = -1. 
      // Check if any of RchSet node has min deg, if so set threshold, SrchPos etc.

      panMarker[nCurrNode] = 0; // Current eliminated supernode representative.
      for(int i = 1; i <= nRchSize; i++)
      {
        const int nRchSetNode = panRchSet[i];
        if(panMarker[nRchSetNode] < 0)
          continue;
        panMarker[nRchSetNode] = 0;
        
        const int nDeg = panDeg[nRchSetNode];
        if(nDeg < nMinDeg) 
          nMinDeg = nDeg;
        if(nDeg <= nThrshld)
        {
          nMinDeg = nThrshld;
          nThrshld = nDeg;
          nSrchPos = panInvPerm[nRchSetNode];
        }
      }

      // Update Quotient Graph.
      if(nNbrhdSize > 0)
        TransformQuotientGraph(
        nCurrNode,
        io_psQGraph,
        &sMinDegData);
    }

  } // End While over nodes.

  // Convert Perm and InvPerm back to C_Style indexing, Indices are shifted back to start from 0.
  // Shift Values.
  for(int i = 1; i <= nNumV; i++)
  {
    o_panPerm[i]--;
  }
  // Shift Pointers.
  o_panPerm++; 
  panInvPerm++;
  delete [] panInvPerm;
}

//==========================
// Static function to find the Reachable set from the node Root
// through the eliminated set S.
// Reach(y,S) = Adj(Nbrhd(y,S) U {y})
// Works only with F_Style indexing.
//==========================
void  REORD_SDSMinDeg::FindReachSet(
  const int i_nVertex,
   REORD_QGraph *io_psQGraph,
   REORD_MinDegData_s * o_psMinDegData)
{
  // Check input.
  assert(i_nVertex > 0 && io_psQGraph && !(io_psQGraph->IsEmpty()) && o_psMinDegData && 
         io_psQGraph->GetStyle() == F_Style && o_psMinDegData->GetStyle() == F_Style);

  int nRoot = i_nVertex;
  int nNumV = io_psQGraph->GetNumV();
  assert(nNumV > 0);
  const int * __restrict panXAdj = io_psQGraph->GetXAdj();
  assert(panXAdj);
  const int * __restrict panAdj = io_psQGraph->GetAdj();
  assert(panAdj);
  const int * __restrict panDeg = o_psMinDegData->GetDeg();
  assert(panDeg);
  int * __restrict panMarker = o_psMinDegData->GetMarker(); // Marker changes to indicate Nodes in Reach set.
  assert(panMarker);

  int * __restrict panRchSet = o_psMinDegData->GetRchSet();
  assert(panRchSet);
  int * __restrict panNbrhdSet = o_psMinDegData->GetNbrhdSet();
  assert(panNbrhdSet);
 
  // nRoot must be an active (live) node i.e. Deg[Root] >= 0, return if already eliminated
  if(panDeg[nRoot] < 0) return; 

  // Initialization of Reachable and Neighborhood sets.
  int  nRchSize = 0;
  int  nNbrhdSize = 0; 


  //panMarker[nRoot] = 1; // Temporarily set Root marker to 1 to avoid including it in Reachset. // TODO: Move it outside. ReConsider in QMupdate where it's set to 2.
  int nRtStart = panXAdj[nRoot];
  int nRtEnd = panXAdj[nRoot + 1];

  for(int k = nRtStart; k < nRtEnd; k++)
  {
    int nNabor = panAdj[k];

    // Connection to Node.
    if(panMarker[nNabor] == 0 && panDeg[nNabor] >= 0)
    {
      // Set Marker to 1 i.e. added to Reach set.
      nRchSize++;
      panMarker[nNabor] = 1;
      panRchSet[nRchSize] = nNabor;
    }

    // Connection to eNode.
    else if(panMarker[nNabor] == 0 && panDeg[nNabor] < 0)
    {
      // Marker -1, later this will be considered as merged under Root supernode.
      nNbrhdSize++;
      panMarker[nNabor] = -1;
      panNbrhdSet[nNbrhdSize] = nNabor;

        // Scan the Adj[Nbrhd].
        int nENodeStart = panXAdj[nNabor];
        int nENodeEnd = panXAdj[nNabor + 1];

        int l = nENodeStart;
        while (l < nENodeEnd)
        {
          int nAdjNbr = panAdj[l];
        
          if(nAdjNbr > 0 && panMarker[nAdjNbr] == 0)
          {
            nRchSize++;
            panMarker[nAdjNbr] = 1;
            panRchSet[nRchSize] = nAdjNbr;
          }

          else if(nAdjNbr < 0)  // Read through linked lists.
          {
            int nLinkNode = -nAdjNbr; 
            l = panXAdj[ nLinkNode];      // l = nENodeStart.
            nENodeEnd = panXAdj[ nLinkNode + 1];
            continue;
          }
          else if (nAdjNbr == 0) // End of List. (F_Style)
          {
            break;
          }

          l++;
        }
    }

  }

  o_psMinDegData->SetRchSizeTo(nRchSize);
  o_psMinDegData->SetNbrhdSizeTo(nNbrhdSize); 

}

//==========================
// Static function to find the Reachable set from the node Root
// through the eliminated set S.
// Reach(y,S) = Adj(Nbrhd(y,S) U {y})
// works with Fortran style 1-indexing
//==========================
void  REORD_SDSMinDeg::FindReachSet(
  const int i_nVertex,
   REORD_QGraph *io_psQGraph,
  int &o_nRchSize,
  int *__restrict o_panRchSet,
  int &o_nNbrhdSize,
  int *__restrict o_panNbrhdSet,
   REORD_MinDegData_s *o_psMinDegData)
{
  // Check input.
  assert(i_nVertex > 0 && io_psQGraph && !(io_psQGraph->IsEmpty()) && o_panRchSet && o_panNbrhdSet && io_psQGraph->GetStyle() == F_Style);
  
  int nRoot = i_nVertex;
  int nNumV = io_psQGraph->GetNumV();
  assert(nNumV > 0);
  const int * __restrict panXAdj = io_psQGraph->GetXAdj();
  assert(panXAdj);
  const int * __restrict panAdj = io_psQGraph->GetAdj();
  assert(panAdj);
  const int * __restrict panDeg = o_psMinDegData->GetDeg();
  assert(panDeg);
  int * __restrict panMarker = o_psMinDegData->GetMarker(); // Marker changes to indicate Nodes in Reach set.
  assert(panMarker);

  // nRoot must be an active (live) node i.e. Deg[Root] >= 0, return if already eliminated
  if(panDeg[nRoot] < 0) return; 

  // Initialization of Reachable and Neighborhood sets.
  o_nRchSize = 0;
  o_nNbrhdSize = 0; 


  //panMarker[nRoot] = 1; // Temporarily set Root marker to 1 to avoid including it in Reachset. // TODO: Move it outside. ReConsider in QMupdate where it's set to 2.
  int nRtStart = panXAdj[nRoot];
  int nRtEnd = panXAdj[nRoot + 1];

  for(int k = nRtStart; k < nRtEnd; k++)
  {
    int nNabor = panAdj[k];

    // Connection to Node.
    if(panMarker[nNabor] == 0 && panDeg[nNabor] >= 0)
    {
      // Set Marker to 1 i.e. added to Reach set.
      o_nRchSize++;
      panMarker[nNabor] = 1;
      o_panRchSet[o_nRchSize] = nNabor;
    }

    // Connection to eNode.
    else if(panMarker[nNabor] == 0 && panDeg[nNabor] < 0)
    {
      // Marker -1, later this will be considered as merged under Root supernode.
      o_nNbrhdSize++;
      panMarker[nNabor] = -1;
      o_panNbrhdSet[o_nNbrhdSize] = nNabor;

        // Scan the Adj[Nbrhd].
        int nENodeStart = panXAdj[nNabor];
        int nENodeEnd = panXAdj[nNabor + 1];

        int l = nENodeStart;
        while (l < nENodeEnd)
        {
          int nAdjNbr = panAdj[l];
        
          if(nAdjNbr > 0 && panMarker[nAdjNbr] == 0)
          {
            o_nRchSize++;
            panMarker[nAdjNbr] = 1;
            o_panRchSet[o_nRchSize] = nAdjNbr;
          }

          else if(nAdjNbr < 0)  // Read through linked lists.
          {
            int nLinkNode = -nAdjNbr; 
            l = panXAdj[ nLinkNode];      // l = nENodeStart.
            nENodeEnd = panXAdj[ nLinkNode + 1];
            continue;
          }
          else if (nAdjNbr == 0) // End of List. (F_Style)
          {
            break;
          }

          l++;
        }
    }

  }

}

//=============================
// static function for the quotient graph transformation after 
// elimination of a node or xnode
// First Reachset is added to the representative of new supernode,
// Then this representative is added to the adj list of Reachset nodes.
//=============================
void  REORD_SDSMinDeg::TransformQuotientGraph(
  const int i_nVertex,
   REORD_QGraph *__restrict io_psQGraph,
  const  REORD_MinDegData_s * __restrict i_psMinDegData)
{
  // Check input.
  assert(i_nVertex > 0 && io_psQGraph && !(io_psQGraph->IsEmpty()) && i_psMinDegData &&
    io_psQGraph->GetStyle() == F_Style && i_psMinDegData->GetStyle() == F_Style);

  int nNumV = io_psQGraph->GetNumV();
  assert(nNumV > 0);
  const int * __restrict panXAdj = io_psQGraph->GetXAdj();
  assert(panXAdj);
  int * __restrict panAdj = io_psQGraph->GetAdj();
  assert(panAdj);
  const int * __restrict panDeg = i_psMinDegData->GetDeg();
  assert(panDeg);
  const int * __restrict panMarker = i_psMinDegData->GetMarker();  // Marker changes to indicate Nodes in Reach set.
  assert(panMarker);

  const int * __restrict panRchSet = i_psMinDegData->GetRchSet();
  assert(panRchSet);
  const int * __restrict panNbrhdSet = i_psMinDegData->GetNbrhdSet();
  assert(panNbrhdSet);
  const int  nRchSize = i_psMinDegData->GetRchSize();
  const int  nNbrhdSize = i_psMinDegData->GetNbrhdSize();


  int nNode = i_nVertex; 
  int nRchCount = 0;
  int nNbrCount = 0;

  // Put ReachSet in the Adj of Node


  bool bRunFlag = true;

  while(true)
  {

    int nStart = panXAdj[nNode];
    int nStop = panXAdj[nNode + 1] - 2; // keep the last one empty for linking.


    if(nStop >= nStart)
    {
      // start inserting ReachSet Nodes into Adj
      for(int i = nStart; i <= nStop; i++)
      {
        nRchCount++;
        panAdj[i] = panRchSet[nRchCount];
    
        if(nRchCount == nRchSize)
        {
          panAdj[i + 1] = 0; // Close the list.
          bRunFlag = false;
          break;
        }
      }
    }
    
    if(!bRunFlag)
      break;

    int nLink = panAdj[nStop + 1];
    nNode = - nLink;
    if(nLink >= 0)
    {
      nNbrCount++;
      nNode = panNbrhdSet[nNbrCount];
      panAdj[nStop + 1] = - nNode;
    }

  }

  // Add the representative eliminated node to the adjacency of nodes in Reach set.

  for(int i = 1; i <= nRchSize; i++)
  {
    const int nRchNode = panRchSet[i];
    const int nStart = panXAdj[nRchNode];
    const int nStop = panXAdj[nRchNode + 1] - 1;

    for(int j = nStart; j <= nStop; j++)
    {
      const int nNabor = panAdj[j];
      if(panMarker[nNabor] < 0)
      {
        panAdj[j] = i_nVertex;      // Place root in place of neighbour nodes.
        break;
      }
    }

  }
   
}


//============================
// static function for updating degrees and merging indistinguishable supernodes.
//
// Works with F_Style.
//============================
void  REORD_SDSMinDeg::UpdateDeg(
     REORD_QGraph * __restrict io_psQGraph,
     REORD_MinDegData_s * i_psMinDegData)
{
  // Check input.
  assert(io_psQGraph && i_psMinDegData && io_psQGraph->GetStyle() == F_Style && i_psMinDegData->GetStyle() == F_Style);
  
  int nNumV = io_psQGraph->GetNumV();
  assert(nNumV > 0);
  int * __restrict panXAdj = io_psQGraph->GetXAdj();
  assert(panXAdj);
  int * __restrict panAdj = io_psQGraph->GetAdj();
  assert(panAdj);
  int * __restrict panDeg = i_psMinDegData->GetDeg();
  assert(panDeg);
  int * __restrict panMarker = i_psMinDegData->GetMarker();  // Marker changes to indicate Nodes in Reach set.
  assert(panMarker);
  int * __restrict panQLink = i_psMinDegData->GetQLink();
  assert(panQLink);
  int * __restrict panQSize = i_psMinDegData->GetQSize();
  assert(panQSize);
  
  int *  panRchSet = i_psMinDegData->GetRchSet();
  assert(panRchSet);
  int *  panNbrhdSet = i_psMinDegData->GetNbrhdSet();
  assert(panNbrhdSet);
  const int  nRchSize = i_psMinDegData->GetRchSize();
  assert(nRchSize > 0);
  const int  nNbrhdSize = i_psMinDegData->GetNbrhdSize();


  int  nDeg0 = 0;               // Cumulatively add the number of the Nodes/Xnodes in RchSet, later to be used for Deg updates.
  int  nNumElimSupNodes = 0;    // Cumulatively add the number of Eliminated Super nodes.

  // If no node in RchSet return.
  if(nRchSize <= 0)
    return;

  // For work arrays use the rest of RchSet and NbrhdSet beyond RchSize and NbrhdSize.
  int * panElimSupNode = &(panNbrhdSet[nNbrhdSize]); // To store eliminated supernodes (C_i's) adjacent to some nodes in Reach Set.
  nDeg0 = 0;
  nNumElimSupNodes = 0;

  //---- Find Eliminated supernodes adjacent to some nodes in ReachSet.
  for(int i = 1; i <= nRchSize; i++)
  {
    const int nNode = panRchSet[i];
    nDeg0 += panQSize[nNode];
    const int nStart = panXAdj[nNode];
    const int nStop = panXAdj[nNode + 1] - 1;

    for(int j = nStart; j <= nStop; j++)
    {
      const int nNabor = panAdj[j];
      if(panMarker[nNabor] == 0 && panDeg[nNabor] < 0)
      {
        panMarker[nNabor] = -1; // Temporarily set to -1 to avoid duplicate.
        nNumElimSupNodes++;
        panElimSupNode[nNumElimSupNodes] = nNabor;
      }
    }
  } // End for over RchSet.

  i_psMinDegData->SetRchTotNodesNumTo(nDeg0);
  i_psMinDegData->SetNumElimSupNodesTo(nNumElimSupNodes);

  //---- Reset Marker for ElimSupNodes and merge if applicable.
  if(nNumElimSupNodes > 0)
  {
    // Reset the Marker value for ElimSupNodes to 0.
    for(int i = 1; i <=nNumElimSupNodes; i++)
      panMarker[panElimSupNode[i]] = 0;   
   
    // Merge XNodes if any.
    MergeXNodes(
      io_psQGraph, 
      i_psMinDegData);
  }

  //---- Update Deg for Nodes in Rch set and not already merged in XNodes.
  
  // Temporary work arrays used to update Deg of simple Nodes in Rch set. 
  // Rchset and Nbrhdset beyond RchSize and NbrhdSize are used for the temporary Rchset and Nbrhdset. 
  int *panTmpNbrhdSet = &(panNbrhdSet[nNbrhdSize]); 
  int *panTmpRchSet = &(panRchSet[nRchSize]);

  for(int i = 1; i <= nRchSize; i++)
  {
    const int nNode = panRchSet[i];
    const int nMark = panMarker[nNode];

    if(nMark == 1) // Maybe nMark == 0, too! 2 and -1 are XNodes already updated.
    {
      panMarker[nNode] = 2; // All Nodes with updated degrees will have Marker = 2.
      int nTmpRchSize = 0;
      int nTmpNbrhdSize = 0;

      FindReachSet(
        nNode,
        io_psQGraph,
        nTmpRchSize,
        panTmpRchSet,
        nTmpNbrhdSize,
        panTmpNbrhdSet,
        i_psMinDegData);

      // Go through tmp Rch set and count their nodes
      int nDeg1 = nDeg0;
      for(int j = 1; j <= nTmpRchSize; j++)
      {
        const int nTmpNode = panTmpRchSet[j];
        nDeg1 += panQSize[nTmpNode];
        panMarker[nTmpNode] = 0;
      }

      // Update Deg.
      panDeg[nNode] = nDeg1 - 1;

      // Reset Tmp Nbrhd markers.
      for(int k = 1; k <= nTmpNbrhdSize; k ++)
      {
        const int nTmpNode = panTmpNbrhdSet[k];
        panMarker[nTmpNode] = 0;
      }

    }
  }
}


//============================
// static function for finding and merging indistinguishable supernodes (XNodes).
// For every y in Y={y in Intersect(R2,R1) | Adj(y) C R1+R2+C1+C2} y's are indistinguishable.
// Works with F_Style.
//============================
void  REORD_SDSMinDeg::MergeXNodes(
   REORD_QGraph * __restrict io_psQGraph,
   REORD_MinDegData_s * __restrict i_psMinDegData)
{
  // Check input.
  assert(io_psQGraph && i_psMinDegData && io_psQGraph->GetStyle() == F_Style && i_psMinDegData->GetStyle() == F_Style);

  int * __restrict panXAdj = io_psQGraph->GetXAdj();
  assert(panXAdj);
  int * __restrict panAdj = io_psQGraph->GetAdj();
  assert(panAdj);
  int * __restrict panDeg = i_psMinDegData->GetDeg();
  assert(panDeg);
  int * __restrict panMarker = i_psMinDegData->GetMarker();  // Marker changes to indicate Nodes in Reach set.
  assert(panMarker);
  int * __restrict panQLink = i_psMinDegData->GetQLink();
  assert(panQLink);
  int * __restrict panQSize = i_psMinDegData->GetQSize();
  assert(panQSize);
  
  int * __restrict panRchSet = i_psMinDegData->GetRchSet();
  assert(panRchSet);
  int * __restrict panNbrhdSet = i_psMinDegData->GetNbrhdSet();
  assert(panNbrhdSet);

  const int nRchSize = i_psMinDegData->GetRchSize();
  const int nNbrhdSize = i_psMinDegData->GetNbrhdSize();
  const int nDeg0 = i_psMinDegData->GetRchTotNodesNum(); // !!
  const int nNumElimSupNodes = i_psMinDegData->GetNumElimSupNodes();

  // Return if there are no Eliminated Supernodes adjacent.
  if(nNumElimSupNodes <= 0)
    return;

  // Work arrays.
  const int * __restrict panElimSupNode = &(panNbrhdSet[nNbrhdSize]); // Nbrhd array is segmented into three parts: (i) Nbrhd, (ii) ElimSupNodes, (iii) Overlap i.e. intersect(R1,R2)
  int * __restrict panOverlap = &(panNbrhdSet[nNbrhdSize + nNumElimSupNodes]); // Overlap for the intersection of Rch set and that of R2=Adj(C_2) .
  int * __restrict panDiffNodes = &(panRchSet[nRchSize]);  // Rch set beyond Rch Size is used for R2-R1.

  
  // Loop over each ElimSupNode, determine R2-R1, Intersect(R2,R1), Check XNode condition, merge if any.
  for(int i = 1; i <= nNumElimSupNodes; i++)
  {
    const int nElimNode = panElimSupNode[i];
    panMarker[nElimNode] = -1;  // To avoid violation when checking the XNode condition will be set to  0 in the end.
    int nNumOvlp = 0;           // Number of nodes in overlap i.e. Intersect(R1,R_i).
    int nNumDiff = 0;           // Number of nodes in R_i - R1.
    int nDeg1= 0;               // Number of nodes in R_i - R1.
    
    int nStart = panXAdj[nElimNode];
    int nStop = panXAdj[nElimNode + 1] - 1;
    int l = nStart;

    // Determine Diff (Ri-R1) and Overlap (Intersect(R_i, R1).
    while(l <= nStop)
    {

      const int nNabor = panAdj[l];

      if(nNabor > 0)
      {
        int nMark = panMarker[nNabor];
        if(nMark == 0) // Diff, R_i -R1, i>1
        {
          nNumDiff++;
          panDiffNodes[nNumDiff] = nNabor;
          nDeg1 += panQSize[nNabor]; 
          panMarker[nNabor] = 1;
        }
        else if (nMark == 1) // If mark < 0 (i.e. some node in previous XNode) it automatically continues until reaches the head of that XNode.
        {
          nNumOvlp++;
          panOverlap[nNumOvlp] = nNabor;
          panMarker[nNabor] = 2;
        }
      }

      else if (nNabor < 0)
      {
        const int nLinkNode = -nNabor;
        l = panXAdj[nLinkNode];
        nStop = panXAdj[nLinkNode + 1] -1;
        continue;
      }

      else if (nNabor == 0) // End of list.
      {
        break;
      }

      l++;

    } // End while Loop over ElimNode Adj.

    // Check the XNode condition and merge if satisfied.
    int nHead = 0;
    int nMrgSize = 0; // Number of nodes in the XNode. 

    for (int i = 1; i <= nNumOvlp; i++)
    {
      const int nOvlpNode = panOverlap[i];
      const int nStart = panXAdj[nOvlpNode];
      const int nStop = panXAdj[nOvlpNode + 1] - 1;
      bool bMrgFlag = true;

      for (int j = nStart; j <= nStop; j++)
      {
        const int nNabor = panAdj[j];
        if (panMarker[nNabor] == 0)    // nNabor not in R1+R2+C1+C2
        {
          panMarker[nOvlpNode] = 1;
          bMrgFlag = false;
          break;
        }
      }

      // Merge if XNode condition satisfied.
      if (bMrgFlag)
      {
        nMrgSize += panQSize[nOvlpNode];
        panMarker[nOvlpNode] = -1;       // Merged. If Head will be set to 2 later.    
        int nNode = nOvlpNode;
        int nLinkNode; 
        while (true)
        {
          nLinkNode = panQLink[nNode];
          if (nLinkNode > 0)
          {
            nNode = nLinkNode;
          }
          else if(nLinkNode == 0)
          {
            panQLink[nNode] = nHead;
            nHead = nOvlpNode;
            break;
          }
        }
      }
    }

    // Update QSize, Deg and Marker
    if(nHead > 0)
    {
      panQSize[nHead] = nMrgSize;
      panDeg[nHead] = nDeg0 + nDeg1 - 1;
      panMarker[nHead] = 2;
    }

    // Reset Marker.
    panMarker[nElimNode] = 0; // Reset Marker of ElimNode from -1 to 0.

    // Reset Nodes in R_i - R1 from 1 to 0.
    for (int j = 1; j <= nNumDiff; j++)
    {
      const int nNode = panDiffNodes[j];
      panMarker[nNode] = 0;
    }

  } // End ElimSupNode loop.

}

//==========================
//
// MMD Multiple Minimum Degree.
//==========================

// Constructor.
  REORD_MMD_s:: REORD_MMD_s(const int i_nNumV)
{
  // Check input.
  assert(i_nNumV > 0);

  try
  {
    m_panHead = new int [i_nNumV];
    m_panQsize = new int [i_nNumV];
    m_panList = new int [i_nNumV];
    m_panMarker = new int [i_nNumV];

    m_nNumV = i_nNumV;
    m_nStyle = C_Style;
  }
  catch(...)
  {
    delete [] m_panHead;
    delete [] m_panQsize;
    delete [] m_panList;
    delete [] m_panMarker;

    assert(false);
  }

}

 // Destructor.
 REORD_MMD_s::~ REORD_MMD_s()
{ 
  printf("\nDeleting MMD\n");

  // Adjusting pointer based on indexing style.
  if(m_nStyle == F_Style)
  {
    m_panHead++;
    m_panQsize++;
    m_panList++;
    m_panMarker++;

  }

  if(m_panHead)
  {
    delete [] m_panHead;
    m_panHead = 0;
  }

  if(m_panQsize)
  {
    delete [] m_panQsize;
    m_panQsize = 0;
  }

  if(m_panList)
  {
    delete [] m_panList;
    m_panList = 0;
  }

  if(m_panMarker)
  {
    delete [] m_panMarker;
    m_panMarker = 0;
  }

  m_panFwd = 0;
  m_panBwd = 0;
}

void  REORD_MMD_s::ConvertToFort()
{
  // Check input.
  assert(!IsEmpty() && GetStyle() == C_Style);
  
  // Shift Pointer.
  m_panHead--;
  m_panQsize--;
  m_panList--;
  m_panMarker--;


  SetStyleTo(F_Style);
}

void  REORD_MMD_s::ConvertToC()
{
  // Check input.
  assert(!IsEmpty() && GetStyle() == F_Style);
  
  // Shift Pointer.
  m_panHead++;
  m_panQsize++;
  m_panList++;
  m_panMarker++;

  SetStyleTo(C_Style);
}

void  REORD_MMD_s::SetFwdBckwdPointers(
  int *i_panInvPerm, 
  int *i_panPerm)
{
  // Check input.
  assert(i_panInvPerm && i_panPerm);

  m_panFwd = i_panInvPerm;
  m_panBwd = i_panPerm;

}

//==============================================================
// Static class for MMD.
// F_Style
// Mehdi P
//==============================================================

//===========================
// MMD Data initialization.
//===========================
void  REORD_SDSMMD::InitDataMMD(
  const  REORD_QGraph *i_psQGraph, 
   REORD_MMD_s *o_psMMDData)
{
  // Check Input.
  assert(i_psQGraph && i_psQGraph->GetStyle() == F_Style &&
         o_psMMDData && o_psMMDData->GetStyle() == F_Style);

  const int nNumV = i_psQGraph->GetNumV();
  const int *  panXAdj= i_psQGraph->GetXAdj(); 
  const int *  panAdj= i_psQGraph->GetAdj();

  int *  panHead = o_psMMDData->GetHead();
  int *  panFwd = o_psMMDData->GetFwd();
  int *  panBwd = o_psMMDData->GetBwd();
  int *  panQsize = o_psMMDData->GetQsize();
  int *  panList = o_psMMDData->GetList();
  int *  panMarker = o_psMMDData->GetMarker();

  // F_Style
  for(int i = 1; i <= nNumV; i++)
  {
    panHead[i] = 0;
    panQsize[i] = 1;
    panList[i] = 0;
    panMarker[i] = 0;
  }

  // Degree doubly linked list.
  for(int nV = 1; nV <= nNumV; nV++)
  {
    const int nDeg = panXAdj[nV + 1] - panXAdj[nV] + 1; // +1 for it count the self edge.
    const int nFV = panHead[nDeg];
    panFwd[nV] = nFV;
    panHead[nDeg] = nV;

    if(nFV > 0)
      panBwd[nFV] = nV;
    panBwd[nV] = -nDeg;
  }
}

//===========================
// Find MMD permutation.
//===========================
void  REORD_SDSMMD::FindMMDPerm(
  const int &i_nDelta,
   REORD_QGraph *io_psQGraph,
  int * o_panPerm,
  int * o_panInvPerm,
  SDS_INT &o_nNumSubs)
{
  // Check input.
  assert(io_psQGraph && !io_psQGraph->IsEmpty() && o_panPerm && o_panInvPerm);

  const int nNumV = io_psQGraph->GetNumV();
  const int nMXINT = INT_MAX;

  // Define MinDegData.
   REORD_MMD_s sMMDData(nNumV); // TODO: check Free Mem.

  // Convert to F_Style.
  io_psQGraph->ConvertToFort();
  sMMDData.ConvertToFort();
  
  o_panPerm--;
  o_panInvPerm--;  
  
  // Use Invperm and Perm for Fwd and Bwd in MMDData. 
  sMMDData.SetFwdBckwdPointers(o_panInvPerm, o_panPerm);


  // Initialize MMDData and fill out the deg list.
  InitDataMMD(
    io_psQGraph, 
    &sMMDData);
  o_nNumSubs = 0;


  int *  panHead = sMMDData.GetHead();
  assert(panHead);
  int *  panMarker = sMMDData.GetMarker();
  assert(panMarker);
  int *  panQsize = sMMDData.GetQsize();
  assert(panQsize);
  int *  panList = sMMDData.GetList();
  assert(panList);

  int nNum = 1;       // Number of ordered nodes.
  int nMDegVtx, nNxtV;  // MinDeg Node (Vertex) and next node in forward.
  int nMDeg = 1;      // MinDeg.
  int nTag;           // Tag is used for marking nodes.

  // Eliminate isolated nodes.
  nNxtV = panHead[nMDeg];
  while(nNxtV > 0)
  {
    nMDegVtx = nNxtV;
    nNxtV = o_panInvPerm[nMDegVtx]; // Forward[nMDegVtx].
    panMarker[nMDegVtx] = nMXINT;
    o_panInvPerm[nMDegVtx] = -nNum;
    nNum++;
  }

  if(nNum > nNumV)
    goto n1000;

  // close Deg structure for isolated nodes, set MDeg for the rest.
  panHead[nMDeg] = 0;
  nMDeg = 2;
  nTag = 1;

  // Loop (1) over all sets of nodes with the same mindeg.
  while(1)
  {
    while(panHead[nMDeg] <=0)
      nMDeg++;

    const int nMDegLmt = nMDeg + i_nDelta;
    int nEHead = 0;

    // Loop (2) within each set of nodes with the same mindeg.
n500:
      nMDegVtx = panHead[nMDeg];
      while(nMDegVtx <= 0)
      {
        nMDeg++;
        if (nMDeg > nMDegLmt)
          goto n900;
        nMDegVtx = panHead[nMDeg];
      }

      // Remove mDegV from the degree struture.
      nNxtV = o_panInvPerm[nMDegVtx];
      panHead[nMDeg] = nNxtV;
      if(nNxtV > 0)
        o_panPerm[nNxtV] = -nMDeg;
      o_panInvPerm[nMDegVtx] = -nNum;
      o_nNumSubs += nMDeg + panQsize[nMDegVtx] -2;
      // Exit if done.
      if((nNum + panQsize[nMDegVtx]) > nNumV)
        goto n1000;

      // Eliminate MinDegNode and Update quotient graph.
      // Reset tag value if necessary.
      nTag++;
      if(nTag >= nMXINT)
      {
        nTag = 1;
        for(int i = 1; i <= nNumV; i++)
        {
          if(panMarker[i] < nMXINT)
            panMarker[i] = 0;
        }
      }

      // Update Quotient graph.
       REORD_SDSMMD::MMDQGraphTransform(
        nMDegVtx,
        io_psQGraph,
        &sMMDData,
        nMXINT,
        nTag);

      nNum += panQsize[nMDegVtx];
      panList[nMDegVtx] = nEHead;
      nEHead = nMDegVtx;

      // if delta < 0, then update happens after every elimination (not a Multiple-MD anymore).
      if(i_nDelta >= 0)
        goto n500;

n900:
      if(nNum >= nNumV)
        goto n1000;

       REORD_SDSMMD::MMDUpdateDeg(
        io_psQGraph,
        &sMMDData,
        nEHead,
        i_nDelta,
        nMDeg,
        nMXINT,
        nTag);

  }

n1000:

 REORD_SDSMMD::MMDNumber(
  nNumV,
  &sMMDData,
  o_panInvPerm,
  o_panPerm);

// Convert to C 0_based indexing.
  io_psQGraph->ConvertToC();
  sMMDData.ConvertToC();

  for(int i = 1; i <= nNumV; i++)
  {
    o_panPerm[i]--;
    o_panInvPerm[i]--;
  }

  o_panInvPerm++;
  o_panPerm++;
}

//===========================
// Eliminate MinDeg node and transform quotient graph.
//===========================
void  REORD_SDSMMD::MMDQGraphTransform(
  const int &i_nMDegVtx,
   REORD_QGraph *io_psQGraph,
   REORD_MMD_s *io_psMMDData,
  const int i_nMXINT,
  const int i_nTag)
{
  // Check input.
  assert(i_nMDegVtx > 0 && io_psQGraph && io_psMMDData && i_nTag < i_nMXINT);

  const int nNumV = io_psQGraph->GetNumV();

  int *  panXAdj = io_psQGraph->GetXAdj();
  assert(panXAdj);
  int *  panAdj = io_psQGraph->GetAdj();
  assert(panAdj);

  int *  panHead = io_psMMDData->GetHead();
  assert(panHead);
  int *  panFwd = io_psMMDData->GetFwd();
  assert(panFwd);
  int *  panBwd = io_psMMDData->GetBwd();
  assert(panBwd);
  int *  panQsize = io_psMMDData->GetQsize();
  assert(panQsize);
  int *  panList = io_psMMDData->GetList();
  assert(panList);
  int *  panMarker = io_psMMDData->GetMarker();
  assert(panMarker);

  // Find Reachable vertexes from mindegnode and place it in the data structure.
  panMarker[i_nMDegVtx] = i_nTag;

  int nStart = panXAdj[i_nMDegVtx];
  int nStop  = panXAdj[i_nMDegVtx + 1] - 1;
  int nLink;

  int nXElimHead = 0;
  int nRchLoc = nStart;
  int nRchLmt = nStop;

  // Keep active nodes of Rch set and form linked list of x-eliminated nodes 
  // (not coincident with nodes of this mindeg set).

  for(int i = nStart; i <= nStop; i++)
  {
    const int nNabor = panAdj[i];
    if(nNabor == 0) 
      break;
    if(panMarker[nNabor] < i_nTag)
    {
      panMarker[nNabor] = i_nTag;  // Tag all the adjacent nodes active and X-eliminated.
      if(panFwd[nNabor] < 0)
      {
        // X-eliminated nodes (neighbors set).
        panList[nNabor] = nXElimHead;
        nXElimHead = nNabor;
      }
      else
      {
        panAdj[nRchLoc] = nNabor;
        nRchLoc++;
      }
    }
  }

  //--- Find Rch nodes via XElimNodes and add them to Adj of current min deg node.
  while(nXElimHead > 0)
  {
    panAdj[nRchLmt] = -nXElimHead;
    nLink = nXElimHead;

n400:
    int nNbrStart = panXAdj[nLink];
    int nNbrStop = panXAdj[nLink + 1] - 1;

    for(int i = nNbrStart; i <= nNbrStop; i++)
    {
      const int nNode = panAdj[i];
      
      nLink = -nNode;
      if(nNode < 0)
        goto n400;
      if(nNode == 0) 
        break;

      if((panMarker[nNode] < i_nTag) && (panFwd[nNode] >= 0))
      {
        panMarker[nNode] = i_nTag;

        if(nRchLoc >= nRchLmt)
        {
          nLink = -panAdj[nRchLmt];
          nRchLoc = panXAdj[nLink];
          nRchLmt = panXAdj[nLink + 1] - 1;
        }
        panAdj[nRchLoc] = nNode;
        nRchLoc++;
      }// End if marker nabor < tag.
    }// End for NbrSrt to NbrStp.

    nXElimHead = panList[nXElimHead];
  }// end of while nXElim > 0.

  // close the new eliminated super node list.
  if(nRchLoc <= nRchLmt)
    panAdj[nRchLoc] = 0;

  //--- For each node in Rch set do the following.
  nLink = i_nMDegVtx;

n1100:
  nStart = panXAdj[nLink];
  nStop = panXAdj[nLink + 1] - 1;

  for(int i = nStart; i <= nStop; i++)
  {
    int nRchNode = panAdj[i];
    nLink = - nRchNode;

    if(nRchNode < 0)
      goto n1100;
    if(nRchNode == 0)
      return;         // no Rch node left return.

    // if not already marked for update or as XNode through elimination of other min deg nodes in the set.
    const int nPvNode = panBwd[nRchNode];
    if((nPvNode != 0) && (nPvNode != -i_nMXINT)) 
    {
      // Deactivate the node in the deg list (either to be updated or in the XNode).
      const int nNxNode = panFwd[nRchNode];

      if(nNxNode > 0) 
        panBwd[nNxNode] = nPvNode;
      if(nPvNode > 0)
        panFwd[nPvNode] = nNxNode;

      const int nDeg = - nPvNode;
      if(nPvNode < 0)
        panHead[nDeg] = nNxNode;
    }// End if Bwd !=0 && Bwd != -Maxint.

    // Distinguish between RchNodes with active (any non-tagged) nodes and else (will become Xnode).
    const int nRchNodeStart = panXAdj[nRchNode];
    const int nRchNodeStop = panXAdj[nRchNode + 1] -1;
    int nActiveNbrCounter = nRchNodeStart; 
    for(int i = nRchNodeStart; i <= nRchNodeStop; i++)
    {
      const int nNabor = panAdj[i];
      if(nNabor == 0)
        break;
      if(panMarker[nNabor] < i_nTag)
      {
        panAdj[nActiveNbrCounter] = nNabor;
        nActiveNbrCounter++;
      }
    }

    // Purge nodes with no active nodes (merge into XNode).
    const int nNumActiveNodes = nActiveNbrCounter - nRchNodeStart;
    if(nNumActiveNodes <= 0) // Xnode
    {
      panQsize[i_nMDegVtx] += panQsize[nRchNode];
      panQsize[nRchNode] = 0;
      panMarker[nRchNode] = i_nMXINT;
      panFwd[nRchNode] = -i_nMDegVtx;
      panBwd[nRchNode] = -i_nMXINT;
    }
    else // Mark to be updated.
    {
      panFwd[nRchNode] = nNumActiveNodes + 1;
      panBwd[nRchNode] = 0;
      panAdj[nActiveNbrCounter] = i_nMDegVtx;
      nActiveNbrCounter++;
      if(nActiveNbrCounter <= nRchNodeStop) 
        panAdj[nActiveNbrCounter] = 0;
    }
  }// End for Strt to Stop.
}

//===========================
// Update degree of updatable nodes.
//===========================
void  REORD_SDSMMD::MMDUpdateDeg(
  const  REORD_QGraph *i_psQGraph,
   REORD_MMD_s *io_psMMDData,
  const int i_nEHead,
  const int i_nDelta,
  int &io_nMinDeg,
  const int i_nMXINT,
  int &io_nTag)
{
  // Check input.
  assert(i_psQGraph && io_psMMDData && i_nEHead > 0 && io_nTag < i_nMXINT);

  const int nNumV = i_psQGraph->GetNumV();

  const int *  panXAdj = i_psQGraph->GetXAdj();
  assert(panXAdj);
  const int *  panAdj = i_psQGraph->GetAdj();
  assert(panAdj);

  int *  panHead = io_psMMDData->GetHead();
  assert(panHead);
  int *  panFwd = io_psMMDData->GetFwd();
  assert(panFwd);
  int *  panBwd = io_psMMDData->GetBwd();
  assert(panBwd);
  int *  panQsize = io_psMMDData->GetQsize();
  assert(panQsize);
  int *  panList = io_psMMDData->GetList();
  assert(panList);
  int *  panMarker = io_psMMDData->GetMarker();
  assert(panMarker);


  const int nMDeg0 = io_nMinDeg + i_nDelta; 
  int nElement = i_nEHead;
  int& nTag = io_nTag;    // Alias to io_nTag.
  int nRNode;

n100:
  // no element left in eliminated node list, return.
  if(nElement <= 0)
    return;

  //-- For each eliminated node perform the following:

  // Reset tag if necessary.
  int nMTag = nTag + nMDeg0;
  if(nMTag >= i_nMXINT)
  {
    nTag = 1;
    for(int i = 1; i <= nNumV; i++)
    {
      if(panMarker[i] < i_nMXINT)
        panMarker[i] = 0;
    }
    nMTag = nTag + nMDeg0;
  }

  //-- Separate Rch nodes into two groups: 
  //   (i) those with two nbrs (Fwd[Rnode] = 2) one active and the other is the eliminated node,
  //   (ii) and those with more than one active nbr i.e., Fwd[Rnode] > 2.
  int nQ2Head = 0;
  int nQXHead = 0;
  int nDeg0 = 0;
  int nLink = nElement;
  int nElimNodeNbr;

n400:
  int iStart = panXAdj[nLink];
  int iStop = panXAdj[nLink + 1] - 1;

  for(int i = iStart; i <= iStop; i++)
  {
    nRNode = panAdj[i];
    nLink = -nRNode;
    if(nRNode < 0)
      goto n400;
    if(nRNode == 0)
      break;

    if(panQsize[nRNode] != 0)
    {
      nDeg0 += panQsize[nRNode];
      panMarker[nRNode] = nMTag;

      // RNode is to be updated (sometimes they get marked as XNode or Outmatched in the previous steps).
      if(panBwd[nRNode] == 0)
      {
        // QX list.
        if(panFwd[nRNode] != 2)
        {
          panList[nRNode] = nQXHead;
          nQXHead = nRNode;
        }
        // Q2 list.
        else
        {
          panList[nRNode] = nQ2Head;
          nQ2Head = nRNode;
        }
      } // End if Bwd == 0.
    } // End if Qsize != 0.
  } // End for istrt istp.

  //-- Deg update and analysis for Q2 nodes
  nRNode = nQ2Head;
  int nQListFlg = 1;

n900:
  if(nRNode <= 0)
    goto n1500;   // end of Q2 list go to QX list.
  if(panBwd[nRNode] != 0)
    goto n2200;  // it might have become XNode or Outmatched, Bwd = -MXINT, get the next Q2 node.

  nTag++;
  int nDeg = nDeg0;

  iStart = panXAdj[nRNode];
  int nNabor = panAdj[iStart];
  if(nNabor == nElement)
    nNabor = panAdj[iStart + 1];
  nLink = nNabor;
  if(panFwd[nNabor] >= 0)
  {
    nDeg += panQsize[nNabor];
    goto n2100;
  }

  // Nabor is eliminated, update deg, or find XNode or Outmatched nodes.
n1000:
  iStart = panXAdj[nLink];
  iStop = panXAdj[nLink + 1] - 1;
  for(int i = iStart; i <= iStop; i++)
  {
    nElimNodeNbr = panAdj[i];
    nLink = -nElimNodeNbr;
    
    if(nElimNodeNbr != nRNode)
    {
      if(nElimNodeNbr < 0)
        goto n1000;
      if(nElimNodeNbr == 0)
        goto n2100;

      if(panQsize[nElimNodeNbr] != 0)
      {
        if(panMarker[nElimNodeNbr] < nTag)
        {
          // Not yet considered.
          panMarker[nElimNodeNbr] = nTag;
          nDeg += panQsize[nElimNodeNbr];
        }
        else
        {
          // Skip already Outmatched nodes (Bwd = - MAXINT). 
          if(panBwd[nElimNodeNbr] == 0) 
          {
            if(panFwd[nElimNodeNbr] == 2)
            {
              // XNode, merge.
              panQsize[nRNode] += panQsize[nElimNodeNbr];
              panQsize[nElimNodeNbr] = 0;
              panMarker[nElimNodeNbr] = i_nMXINT;
              panFwd[nElimNodeNbr] = -nRNode;
              panBwd[nElimNodeNbr] = -i_nMXINT;
            }
            else
            {
              // nElimNodeNbr is outmatched ny nRNode.
              panBwd[nElimNodeNbr] = -i_nMXINT;
            }
          }// End if panBwd[nElimNodeNbr] == 0.
        }// end else.
      } // End if panQsize[nElimNodeNbr] != 0.
    }// End if nElimNodeNbr != nRNode.
  }// End for iStrt iStp for elminated nbrs for Q2 RNodes.

  goto n2100;

n1500:
  //-- Deg update and analysis for QX nodes
  nRNode = nQXHead;
  nQListFlg = 0;

n1600:
  if(nRNode <= 0)
    goto n2300;
  if(panBwd[nRNode] != 0)
    goto n2200; // get the next QX node.
  nTag++;
  nDeg = nDeg0;
  // For each unmarked nabor of RNode, do the following.
  iStart = panXAdj[nRNode];
  iStop = panXAdj[nRNode + 1] - 1;

  for(int i = iStart; i <= iStop; i++)
  {
    nNabor = panAdj[i];
    if (nNabor == 0)
      break;
    if(panMarker[nNabor] < nTag)
    {
      panMarker[nNabor] = nTag;
      nLink = nNabor;
      // if not eliminated update deg.
      if(panFwd[nNabor] >= 0)
        nDeg += panQsize[nNabor];
      else
      {
n1700:
        int jStart = panXAdj[nLink];
        int jStop = panXAdj[nLink + 1] - 1;
        for(int j = jStart; j <= jStop; j++)
        {
          nElimNodeNbr = panAdj[j];
          nLink = -nElimNodeNbr;
          if(nElimNodeNbr < 0)
            goto n1700;
          if(nElimNodeNbr == 0)
            break;
          if(panMarker[nElimNodeNbr] < nTag)
          {
            panMarker[nElimNodeNbr] = nTag;
            nDeg += panQsize[nElimNodeNbr];
          }
        }// End if jStart jStp.
      } // End else.
    }// End if(panMarker[nNabor] < nTag).
  } // End fo iStrt iStp.

n2100:
  // Find external degree and update.
  nDeg = nDeg - panQsize[nRNode] + 1;
  int nFNode = panHead[nDeg];
  panHead[nDeg] = nRNode;
  panFwd[nRNode] = nFNode;
  panBwd[nRNode] = - nDeg;
  if(nFNode > 0)
    panBwd[nFNode] = nRNode;
  if(nDeg < io_nMinDeg)
    io_nMinDeg = nDeg;

n2200:
  // Get the next Rchnode for either Q2 or QX list.
  nRNode = panList[nRNode];
  if(nQListFlg == 1) 
    goto n900;
  goto n1600; // Go to Qx list.

n2300:
  // Get the next Eliminated node.
  nTag = nMTag;
  nElement = panList[nElement];
  goto n100;

}
  
//===========================
// Clculate the perm and inv perm 
// numbering from info already contained in them.
//===========================
void  REORD_SDSMMD::MMDNumber(
  const int i_nNumV,
  const  REORD_MMD_s *i_psMMDData,
  int * io_panInvPerm,
  int * io_panPerm)
{
  //Check input.
  assert(io_panInvPerm && io_panPerm && i_psMMDData);

  const int *  panQsize = i_psMMDData->GetQsize();
  assert(panQsize);

  // Use Perm array as work array and initialize as follows.
  for(int i = 1; i <= i_nNumV; i++)
  {
    const int nQsize = panQsize[i];
    if(nQsize > 0)
      io_panPerm[i] = -io_panInvPerm[i];
    // Merged XNodes.
    if(nQsize <= 0)
      io_panPerm[i] = io_panInvPerm[i];
  }

  for(int nNode = 1; nNode <= i_nNumV; nNode++)
  {
    // if XNode
    if(io_panPerm[nNode] <= 0)
    {
      int nFather = nNode;

      // Find Root.
      while(io_panPerm[nFather] <=0)
        nFather = -io_panPerm[nFather];

      const int nRoot = nFather;
      const int nNum = io_panPerm[nRoot] + 1;
      io_panInvPerm[nNode] = -nNum;
      io_panPerm[nRoot] = nNum;

      // Shorten the XNode tree.
      nFather = nNode;
      int nNxtFthr = -io_panPerm[nFather];
      while(nNxtFthr > 0)
      {
        io_panPerm[nFather] = -nRoot;
        nFather = nNxtFthr;
        nNxtFthr = -io_panPerm[nFather];
      }
    }
  }

  //-- Generate the numbering using info in InvPerm.
  for(int nNode = 1; nNode <= i_nNumV; nNode++)
  {
    const int nNum = -io_panInvPerm[nNode];
    io_panInvPerm[nNode] = nNum;
    io_panPerm[nNum] = nNode;
  }
}

//===========================
// Symbolic factorization. 
// 
//===========================
void  REORD_Solve::DoSymbFact(
  const  REORD_QGraph *i_psGraph,
  const int *i_panPerm,
  const int *i_panInvPerm,
  int *o_panXLnz,
  int *o_panXNzsub,
  int *o_panNzSub,
  int &o_nNumNNZ,
  int &o_nNumSubs, // io in case not enough memory.
  int &o_eError)
{
  // Check input.
  assert(i_psGraph && i_panInvPerm && i_panPerm && o_panNzSub && o_panXNzsub && o_panXLnz);

  const int nNumV = i_psGraph->GetNumV();
  const int *__restrict panXAdj = i_psGraph->GetXAdj();
  assert(panXAdj);
  const int *__restrict panAdj = i_psGraph->GetAdj();
  assert(panAdj);

  // working arrays memory allocation.
  int *__restrict panRchLnk = NULL;
  int *__restrict panMrgLnk = NULL;
  int *__restrict panMarker = NULL;

  try
  {
    panRchLnk = new int [nNumV];
    panMrgLnk = new int [nNumV];
    panMarker = new int [nNumV];
  }
  catch(...)
  {
    delete[] panRchLnk;
    delete[] panMrgLnk;
    delete[] panMarker;
    assert(false);
    o_eError = 11;
    return;
  }

  // Initialization.
  o_eError = 0;
  int nzbeg = 0;       // begining and end position from the previous column.
  int nzend = -1;      // With no compression, current column starts at nzend of last column.
  o_panXLnz[0] = 0;
  int nTag = nNumV;
  int nRchNodePrev, nRchNodeNxt; // To merge R_i in order.

  for(int i = 0; i < nNumV; i++)
  {
    panMrgLnk[i] = -1;
    panMarker[i] = -1;
  }

  // For every column associated with a vertex (node || equation) perform the following:
  // KNZ counts the number of nonzeros in that column. panRchLnk accumulates the row subscripts.

  for(int K = 0; K < nNumV; K++)
  {
    int KNZ = 0;
    int nMrgNode = panMrgLnk[K];
    int nMarkFlg = 0;

    // for mass symbfact.
    panMarker[K] = K;
    if(nMrgNode >= 0)
      panMarker[K] = panMarker[nMrgNode];

    o_panXNzsub[K] = nzend; 

    const int nNode = i_panPerm[K];
    const int nJStart = panXAdj[nNode];
    const int nJStop = panXAdj[nNode + 1] -1;
    
    if(nJStart > nJStop) 
    {
      o_panXLnz[K + 1] = o_panXLnz[K];
      continue;   // when the column is empty (isolated node). !!!!!!!
    }

    //-- Add the structure of A[*,K] below diagonal to RchLnk.
    panRchLnk[K] = nTag; // tag the end of the list.
    for(int j = nJStart; j <= nJStop; j++)
    {
      int nNabor = panAdj[j];
      nNabor = i_panInvPerm[nNabor];
  
      // Adj(x_k) - S_k.
      if(nNabor <= K)
        continue;

      nRchNodePrev = K;
      nRchNodeNxt = panRchLnk[K];
      
      // Nodes are put in RchLnk in order.
      while(nRchNodeNxt <= nNabor)
      {
        nRchNodePrev = nRchNodeNxt;
        nRchNodeNxt = panRchLnk[nRchNodeNxt];
      }

      KNZ++;
      panRchLnk[nRchNodePrev] = nNabor;
      panRchLnk[nNabor] = nRchNodeNxt;
      if(panMarker[nNabor] != panMarker[K])
        nMarkFlg = 1;
    } // End for j.

    //-- Mass Symbolic Elimination.
    int nLmax = 0;
    if(nMarkFlg == 0 && nMrgNode >= 0 && panMrgLnk[nMrgNode] < 0)
    {
      o_panXNzsub[K] = o_panXNzsub[nMrgNode] + 1;
      KNZ = o_panXLnz[nMrgNode + 1] - (o_panXLnz[nMrgNode] + 1);
      goto n1400;
    }

    //-- Add Rchset from MargLnk to RchLnk.
    nMrgNode = panMrgLnk[K];
    while(nMrgNode >= 0)
    {
      const int nMrgNodeNZ = o_panXLnz[nMrgNode + 1] - (o_panXLnz[nMrgNode] + 1); 
      const int nIStart = o_panXNzsub[nMrgNode] + 1;
      const int nIStop = o_panXNzsub[nMrgNode] + nMrgNodeNZ;

      // In case the RchSet(nMrgNode) with Lmax is a superset of all other MrgNodes -> RchSet(K) = RchSet(MergNode of Lmax).
      if(nMrgNodeNZ > nLmax)
      {
        nLmax = nMrgNodeNZ;
        o_panXNzsub[K] = nIStart;
      }

      int nRchM = K;
      for(int i = nIStart; i <= nIStop; i++)
      {
        const int nNabor = o_panNzSub[i];

        nRchNodePrev = nRchM;
        nRchNodeNxt = panRchLnk[nRchM];

        while(nNabor > nRchNodeNxt)
        {
          nRchNodePrev = nRchNodeNxt;
          nRchNodeNxt = panRchLnk[nRchNodeNxt];
        }
        if(nNabor == nRchNodeNxt)
          continue;

        KNZ++;
        panRchLnk[nRchNodePrev] = nNabor;
        panRchLnk[nNabor] = nRchNodeNxt;
        nRchM = nNabor;
      }

      nMrgNode = panMrgLnk[nMrgNode];

    } // End while MrgNode.
    
    //-- if MrgNode with Lmax is the superset of Adj(x_k) and RchSet(other MrgNodes).
    if(KNZ == nLmax)
      goto n1400;

    //-- else for cases of XNodes, or RchSet(K) C RchSet(K-1) or RchSet(K-1) C RchSet(K) do
    int nXNodeSubPos;
    if(nzbeg < nzend)
    {
      nRchNodeNxt = panRchLnk[K];

      for(nXNodeSubPos = nzbeg; nXNodeSubPos <= nzend; nXNodeSubPos++)
      {
        const int nTmp = o_panNzSub[nXNodeSubPos] - nRchNodeNxt;
        if(nTmp < 0)
          continue;
        // matching node found in both RchSet(K) and RchSet(mergnode).
        if(nTmp == 0)
          goto n1000;
        // no match found proceed to copying form RchLnk to nzsub.
        else if(nTmp > 0)
          goto n1200;
      }
      goto n1200;
    }
    goto n1200;

n1000:
    o_panXNzsub[K] = nXNodeSubPos;
    for(int j = nXNodeSubPos; j <= nzend; j++)
    {
      if(o_panNzSub[j] != nRchNodeNxt)
        goto n1200;
      nRchNodeNxt = panRchLnk[nRchNodeNxt];
      // RchSet(K) C RchSet(MrgNode) -> no need to copy, compress.
      if(nRchNodeNxt == nTag)
        goto n1400;
    }
    // RchSet(MrgNode) C RchSet(K) partially compress.
    nzend = nXNodeSubPos - 1;

n1200:
    // Copy structure of L[*,K] into NZSub and XNZSub.
    nzbeg = nzend + 1;
    nzend = nzend + KNZ;
    nRchNodeNxt = K;

    if(nzend > o_nNumSubs)
    {
      o_panNzSub = (int*) realloc(o_panNzSub, 2 * o_nNumSubs * sizeof(*o_panNzSub));
      o_nNumSubs *= 2;
      if(!o_panNzSub)
      {
        o_eError = -11;
        break;
      }
    }

    for(int i = nzbeg; i <= nzend; i++)
    {
      nRchNodeNxt = panRchLnk[nRchNodeNxt];
      o_panNzSub[i] = nRchNodeNxt;
      panMarker[nRchNodeNxt] = K;
    }

    o_panXNzsub[K] = nzbeg;
    panMarker[K] = K;

n1400:

    //  update the vector mrglnk. note column l(*,k) just found
    // is required to determine column l(*,j), where
    //  l(j,k) is the first nonzero in l(*,k) below diagonal.
    
    o_panXLnz[K + 1] = o_panXLnz[K] + KNZ;
    if(KNZ > 1)
    {
      const int nTmpind = o_panXNzsub[K];
      const int nMinNode = o_panNzSub[nTmpind];
      panMrgLnk[K] = panMrgLnk[nMinNode];
      panMrgLnk[nMinNode] = K;
    }
  } // End for K.

  o_nNumNNZ = o_panXLnz[nNumV]; // off-diagonals are stored in compressed form and the last column col[nNumV - 1] is used for the number of nnz.
  o_panXNzsub[nNumV] = o_panXNzsub[nNumV - 1]; // the last position (i.e. nNumV) in both XLnz and XNzSub just repeat the previous values.
  o_nNumSubs = o_panXNzsub[nNumV];

  o_panNzSub = (int*) realloc(o_panNzSub, o_nNumSubs * sizeof(*o_panNzSub));

  delete [] panRchLnk;
  delete [] panMrgLnk;
  delete [] panMarker;

}


void  REORD_Solve::DoSymbFact2(
  const  REORD_QGraph *i_psGraph,
  const int *i_panPerm,
  const int *i_panInvPerm,
  SDS_INT *o_panXLnz,
  SDS_INT *o_panXNzsub,
  SDS_INT **o_panNzSub,
  SDS_INT &o_nNumNNZ,
  SDS_INT &o_nNumSubs, // io in case not enough memory.
  int &o_eError)
{
  // Check input.
  assert(i_psGraph && i_panInvPerm && i_panPerm && o_panNzSub && o_panXNzsub && o_panXLnz && o_nNumSubs > 0);

  const SDS_INT nNumV = i_psGraph->GetNumV();
  const int *__restrict panXAdj = i_psGraph->GetXAdj();
  assert(panXAdj);
  const int *__restrict panAdj = i_psGraph->GetAdj();
  assert(panAdj);

  // working arrays memory allocation.
  SDS_INT *__restrict panRchLnk = NULL;
  SDS_INT *__restrict panMrgLnk = NULL;
  SDS_INT *__restrict panMarker = NULL;

  try
  {
    panRchLnk = new SDS_INT [nNumV];
    panMrgLnk = new SDS_INT [nNumV];
    panMarker = new SDS_INT [nNumV];
  }
  catch(...)
  {
    delete[] panRchLnk;
    delete[] panMrgLnk;
    delete[] panMarker;
    assert(false);
    o_eError = 11;
    return;
  }

  // Initialization.
  o_eError = 0;
  SDS_INT nzbeg = 0;       // begining and end position from the previous column.
  SDS_INT nzend = -1;      // With no compression, current column starts at nzend of last column.
  o_panXLnz[0] = 0;
  SDS_INT nTag = nNumV;
  SDS_INT nRchNodePrev, nRchNodeNxt; // To merge R_i in order.

  for(SDS_INT i = 0; i < nNumV; i++)
  {
    panMrgLnk[i] = -1;
    panMarker[i] = -1;
  }

  // For every column associated with a vertex (node || equation) perform the following:
  // KNZ counts the number of nonzeros in that column. panRchLnk accumulates the row subscripts.

  for(int K = 0; K < nNumV; K++)
  {
    SDS_INT KNZ = 0;
    SDS_INT nMrgNode = panMrgLnk[K];
    SDS_INT nMarkFlg = 0;

    // for mass symbfact.
    panMarker[K] = K;
    if(nMrgNode >= 0)
      panMarker[K] = panMarker[nMrgNode];

    o_panXNzsub[K] = nzend; 

    const SDS_INT nNode = i_panPerm[K];
    const SDS_INT nJStart = panXAdj[nNode];
    const SDS_INT nJStop = panXAdj[nNode + 1] -1;
    
    if(nJStart > nJStop) 
    {
      o_panXLnz[K + 1] = o_panXLnz[K];
      continue;   // when the column is empty (isolated node). !!!!!!!
    }

    //-- Add the structure of A[*,K] below diagonal to RchLnk.
    panRchLnk[K] = nTag; // tag the end of the list.
    for(SDS_INT j = nJStart; j <= nJStop; j++)
    {
      int nNabor = panAdj[j];
      nNabor = i_panInvPerm[nNabor];
  
      // Adj(x_k) - S_k.
      if(nNabor <= K)
        continue;

      nRchNodePrev = K;
      nRchNodeNxt = panRchLnk[K];
      
      // Nodes are put in RchLnk in order.
      while(nRchNodeNxt <= nNabor)
      {
        nRchNodePrev = nRchNodeNxt;
        nRchNodeNxt = panRchLnk[nRchNodeNxt];
      }

      KNZ++;
      panRchLnk[nRchNodePrev] = nNabor;
      panRchLnk[nNabor] = nRchNodeNxt;
      if(panMarker[nNabor] != panMarker[K])
        nMarkFlg = 1;
    } // End for j.

    //-- Mass Symbolic Elimination.
    SDS_INT nLmax = 0;
    if(nMarkFlg == 0 && nMrgNode >= 0 && panMrgLnk[nMrgNode] < 0)
    {
      o_panXNzsub[K] = o_panXNzsub[nMrgNode] + 1;
      KNZ = o_panXLnz[nMrgNode + 1] - (o_panXLnz[nMrgNode] + 1);
      goto n1400;
    }

    //-- Add Rchset from MargLnk to RchLnk.
    nMrgNode = panMrgLnk[K];
    while(nMrgNode >= 0)
    {
      const SDS_INT nMrgNodeNZ = o_panXLnz[nMrgNode + 1] - (o_panXLnz[nMrgNode] + 1); 
      const SDS_INT nIStart = o_panXNzsub[nMrgNode] + 1;
      const SDS_INT nIStop = o_panXNzsub[nMrgNode] + nMrgNodeNZ;

      // In case the RchSet(nMrgNode) with Lmax is a superset of all other MrgNodes -> RchSet(K) = RchSet(MergNode of Lmax).
      if(nMrgNodeNZ > nLmax)
      {
        nLmax = nMrgNodeNZ;
        o_panXNzsub[K] = nIStart;
      }

      SDS_INT nRchM = K;
      for(SDS_INT i = nIStart; i <= nIStop; i++)
      {
        const SDS_INT nNabor = (*o_panNzSub)[i];

        nRchNodePrev = nRchM;
        nRchNodeNxt = panRchLnk[nRchM];

        while(nNabor > nRchNodeNxt)
        {
          nRchNodePrev = nRchNodeNxt;
          nRchNodeNxt = panRchLnk[nRchNodeNxt];
        }
        if(nNabor == nRchNodeNxt)
          continue;

        KNZ++;
        panRchLnk[nRchNodePrev] = nNabor;
        panRchLnk[nNabor] = nRchNodeNxt;
        nRchM = nNabor;
      }

      nMrgNode = panMrgLnk[nMrgNode];

    } // End while MrgNode.
    
    //-- if MrgNode with Lmax is the superset of Adj(x_k) and RchSet(other MrgNodes).
    if(KNZ == nLmax)
      goto n1400;

    //-- else for cases of XNodes, or RchSet(K) C RchSet(K-1) or RchSet(K-1) C RchSet(K) do
    SDS_INT nXNodeSubPos;
    if(nzbeg < nzend)
    {
      nRchNodeNxt = panRchLnk[K];

      for(nXNodeSubPos = nzbeg; nXNodeSubPos <= nzend; nXNodeSubPos++)
      {
        const SDS_INT nTmp = (*o_panNzSub)[nXNodeSubPos] - nRchNodeNxt;
        if(nTmp < 0)
          continue;
        // matching node found in both RchSet(K) and RchSet(mergnode).
        if(nTmp == 0)
          goto n1000;
        // no match found proceed to copying form RchLnk to nzsub.
        else if(nTmp > 0)
          goto n1200;
      }
      goto n1200;
    }
    goto n1200;

n1000:
    o_panXNzsub[K] = nXNodeSubPos;
    for(SDS_INT j = nXNodeSubPos; j <= nzend; j++)
    {
      if((*o_panNzSub)[j] != nRchNodeNxt)
        goto n1200;
      nRchNodeNxt = panRchLnk[nRchNodeNxt];
      // RchSet(K) C RchSet(MrgNode) -> no need to copy, compress.
      if(nRchNodeNxt == nTag)
        goto n1400;
    }
    // RchSet(MrgNode) C RchSet(K) partially compress.
    nzend = nXNodeSubPos - 1;

n1200:
    // Copy structure of L[*,K] into NZSub and XNZSub.
    nzbeg = nzend + 1;
    nzend = nzend + KNZ;
    nRchNodeNxt = K;

    if(nzend > o_nNumSubs)
    {
      (*o_panNzSub) = (SDS_INT*) realloc((*o_panNzSub), 2 * o_nNumSubs * sizeof(**o_panNzSub));
      o_nNumSubs *= 2;
      if(!(*o_panNzSub))
      {
        o_eError = -11;
        break;
      }
    }

    for(SDS_INT i = nzbeg; i <= nzend; i++)
    {
      nRchNodeNxt = panRchLnk[nRchNodeNxt];
      (*o_panNzSub)[i] = nRchNodeNxt;
      panMarker[nRchNodeNxt] = K;
    }

    o_panXNzsub[K] = nzbeg;
    panMarker[K] = K;

n1400:

    //  update the vector mrglnk. note column l(*,k) just found
    // is required to determine column l(*,j), where
    //  l(j,k) is the first nonzero in l(*,k) below diagonal.
    
    o_panXLnz[K + 1] = o_panXLnz[K] + KNZ;
    if(KNZ > 1)
    {
      const SDS_INT nTmpind = o_panXNzsub[K];
      const SDS_INT nMinNode = (*o_panNzSub)[nTmpind];
      panMrgLnk[K] = panMrgLnk[nMinNode];
      panMrgLnk[nMinNode] = K;
    }
  } // End for K.

  o_nNumNNZ = o_panXLnz[nNumV]; // off-diagonals are stored in compressed form and the last column col[nNumV - 1] is used for the number of nnz.
  o_panXNzsub[nNumV] = o_panXNzsub[nNumV - 1]; // the last position (i.e. nNumV) in both XLnz and XNzSub just repeat the previous values.
  o_nNumSubs = o_panXNzsub[nNumV];

  *o_panNzSub = (SDS_INT*) realloc((*o_panNzSub), o_nNumSubs * sizeof(**o_panNzSub));

  delete [] panRchLnk;
  delete [] panMrgLnk;
  delete [] panMarker;

}

//==============================================================
// Static class for ND.
// Mehdi P
//==============================================================

//===========================
// Level structure..
//===========================
void  REORD_NDOrdering::FindLevelStruct(
  const int &i_nRoot,
  const int *__restrict i_panXAdj,
  const int *__restrict i_panAdj,
  int *__restrict io_panMask,
  int &o_nNumLvl,
  int *__restrict o_panXLS,
  int *__restrict o_panLS)
{
  //Check input.
  assert(i_nRoot >= 0 && i_panXAdj && i_panAdj && io_panMask && o_panXLS && o_panLS);

  // Initialization.
  o_nNumLvl = 0;
  int nLvlStrt;
  int nLvlStp = -1;
  int nLvlNum = -1;
  int nLvlWidth;

  int nCounter = 0;
  o_panLS[nCounter] = i_nRoot;
  io_panMask[i_nRoot] = -1;

  do
  {  
    nLvlStrt = nLvlStp + 1;
    nLvlStp = nCounter;
    nLvlNum++;
    o_panXLS[nLvlNum] = nLvlStrt;  // XLS[n+1] 1 extra space for nnz (here number of nodes in connected comp.

    for(int i = nLvlStrt; i <= nLvlStp; i++)  //!!!!! <=
    {
      const int nNode = o_panLS[i];
      const int nStart = i_panXAdj[nNode];
      const int nStop = i_panXAdj[nNode + 1] -1;
      if(nStop < nStart)
        break; //!!!!!!

      for(int j = nStart; j <= nStop; j++)
      {
        const int nNbr = i_panAdj[j];
        if(io_panMask[nNbr] > -1)
        {
          nCounter++;
          o_panLS[nCounter] = nNbr;
          io_panMask[nNbr] = -1;
        }
      }
    } 

    nLvlWidth = nCounter - nLvlStp;
 
  } while(nLvlWidth > 0);

  o_nNumLvl = nLvlNum + 1;     // note that level length = nLvlNum {L0,L1,..., Le} where level length = e.
  o_panXLS[o_nNumLvl] = nCounter + 1;

  // Reset mask.
  for(int i = 0; i <= nCounter; i++)
  {
    const int nNode = o_panLS[i];
    io_panMask[nNode] = 0;
  }
}

//===========================
// Find Pseudo-peripheral node.
// The section subgraph is determined
// via Root and Mask.
//===========================
void  REORD_NDOrdering::PsedoPeripheralNodeLS(
  int &io_nRoot,
  const int *__restrict i_panXAdj,
  const int *__restrict i_panAdj,
  int *__restrict i_panMask,
  int &o_nNumLvl,
  int *__restrict o_panXLS,
  int *__restrict o_panLS)
{
  //Check input.
  assert(io_nRoot >= 0 && i_panXAdj && i_panAdj  && i_panMask && o_panXLS && o_panLS);

  // Start from an arbitrary node for the starting point (i.e. io_nRoot).
  // Determine the level structure.

  FindLevelStruct(
    io_nRoot,
    i_panXAdj,
    i_panAdj,
    i_panMask,
    o_nNumLvl,
    o_panXLS,
    o_panLS);

   int nNumLvlTotNodes = o_panXLS[o_nNumLvl];

   // if graph is a chain or isolated node, return.
   if(o_nNumLvl == nNumLvlTotNodes || o_nNumLvl == 1)
     return;

   // Find node of mindeg in the last level and establish level structure.
   int nPrevNumLvl;
   do
   {
     nPrevNumLvl = o_nNumLvl;
     const int nJStart = o_panXLS[o_nNumLvl - 1];
     const int nJStop =  nNumLvlTotNodes - 1;
     int nMinDeg = nNumLvlTotNodes;
     io_nRoot = o_panLS[nJStart];

     if(nJStop > nJStart)
     {
       for(int j = nJStart; j <= nJStop; j++)
       {
         const int nNode = o_panLS[j];
         const int nKStrat = i_panXAdj[nNode];
         const int nKStop = i_panXAdj[nNode + 1] -1;
         int nDeg = 0;

         for(int k = nKStrat; k <= nKStop; k++)
         {
           const int nNbr = i_panAdj[k];
           if(i_panMask[nNbr] > -1)
             nDeg++;
         }

         if(nDeg < nMinDeg)
         {
           nMinDeg = nDeg;
           io_nRoot = nNode;
         }
       }
     }

    FindLevelStruct(
      io_nRoot,
      i_panXAdj,
      i_panAdj,
      i_panMask,
      o_nNumLvl,
      o_panXLS,
      o_panLS);

    if(o_nNumLvl == nNumLvlTotNodes)
      return;

   } while(o_nNumLvl > nPrevNumLvl);

}

//===========================
// Nested Diesection permutation.
//===========================
//subroutine gennd ( neqns, xadj, adjncy, mask,
//     &  perm, xls, ls )

void  REORD_NDOrdering::DoNDOrdering(
   REORD_QGraph *i_psGraph,
  int *o_panPerm,
  int *o_panInvPerm)
{
  //Check input.
  assert(i_psGraph && o_panPerm && o_panInvPerm);

  const int nNumV = i_psGraph->GetNumV();
  int *__restrict panXAdj = i_psGraph->GetXAdj();
  assert(panXAdj);
  const int *__restrict panAdj = i_psGraph->GetAdj();
  assert(panAdj);

  // Allocations.
  int *__restrict panMask = NULL;
  int *__restrict panXLS = NULL;
  int *__restrict panLS = NULL;

  try
  {
    panMask = new int [nNumV];
    panXLS = new int [nNumV + 1];
    panLS = new int [nNumV];
  }
  catch(...)
  {
    delete[] panMask;
    delete[] panXLS;
    delete[] panLS;
    assert(false);
    throw;
  }

  // Initialization.
  ::memset(panMask, 0, sizeof(*panMask) * nNumV);
  int nNum = 0;

  // Loop over all nodes and find separator for unmasked connected componentes.
  for(int i = 0; i < nNumV; i++)
  {
    int nRoot;
    if(panMask[i] < 0)
      continue;

    nRoot = i;
    int nSepSize;

    // find and number separator.
    FindSeparator(
      nRoot,
      panXAdj,
      panAdj,
      panMask,
      nSepSize,
      (o_panPerm + nNum),
      panXLS,
      panLS);

    nNum += nSepSize;

    if(nNum >= nNumV)
      break;
  }

  // Reverese order: first separator found should be last.

  ReverseArray(nNumV, o_panPerm);

  // create inv perm array.
  for(int i = 0; i < nNumV; i++)
  {
    const int nIndxOld = o_panPerm[i];
    o_panInvPerm[nIndxOld] = i;
  }
}

//===========================
// Find Separator and number.
//===========================

void  REORD_NDOrdering::FindSeparator(
  int &i_nRoot,
  int *__restrict i_panXAdj,
  const int *__restrict i_panAdj,
  int *__restrict io_panMask,
  int &o_nSepSize,
  int *__restrict o_panSepVec,
  int *__restrict panXLS,
  int *__restrict panLS)
{
  // Check input.
  assert(i_nRoot >= 0 && i_panXAdj && i_panAdj && io_panMask && o_panSepVec && panXLS && panLS);

  o_nSepSize = 0;
  int nNumLvl = 0;
  int nNode;

  // Find level structure rooted at Pseudo-peripheral node.
  PsedoPeripheralNodeLS(
    i_nRoot,
    i_panXAdj,
    i_panAdj,
    io_panMask,
    nNumLvl,
    panXLS,
    panLS);

  // If the number of levels is less than 4 (i.e., L0,L1,L2) , return the whole component as the separator.
  if(nNumLvl <= 3)
  {
    o_nSepSize = panXLS[nNumLvl];
    for(int i = 0; i < o_nSepSize; i++)
    {
      nNode = panLS[i];
      o_panSepVec[i] = nNode;
      io_panMask[nNode] = -1;
    }
    return;
  }

  // Find middle level Lj, form separator by including nodes with connection to the next level Lj+1.

  const int nMidLvlIndx = (nNumLvl - 1) / 2;
  const int nMidLStart = panXLS[nMidLvlIndx];
  const int nNxtLStart = panXLS[nMidLvlIndx + 1];
  const int nMidLStop = nNxtLStart - 1;
  const int nNxtLStop = panXLS[nMidLvlIndx + 2] - 1;

 // The separator is obtained by including only those
 // MIDDLE-level nodes with neighbors in the MIDDLE+1
 // level.  XADJ is used temporarily to mark those
 // nodes in the MIDDLE+1 level.

  for(int i = nNxtLStart; i <= nNxtLStop; i++)
  {
    nNode = panLS[i];
    i_panXAdj[nNode] = - i_panXAdj[nNode];
  }

  for(int i = nMidLStart; i <= nMidLStop; i++)
  {
    nNode = panLS[i];
    const int nStrt = i_panXAdj[nNode];
    const int nStp = abs(i_panXAdj[nNode + 1]) - 1;

    for(int j = nStrt; j <= nStp; j++)
    {
      const int nNbr = i_panAdj[j];
      if(i_panXAdj[nNbr] < 0)
      {
        o_panSepVec[o_nSepSize] = nNode;
        io_panMask[nNode] = -1;
        o_nSepSize++;
        break;
      }
    }
  }

  // return the panXadj to its original form.
  for(int i = nNxtLStart; i <= nNxtLStop; i++)
  {
    nNode = panLS[i];
    i_panXAdj[nNode] = - i_panXAdj[nNode];
  }
}

//===========================
// Reverse integer array.
//===========================
void  REORD_NDOrdering::ReverseArray(
  const int i_nSize,
  int *io_panArray)
{
  // Check input.
  assert(i_nSize > 0 && io_panArray);

  for(int i = 0; i < (i_nSize/2); i++)
  {
    const int nTmp = io_panArray[i];
    io_panArray[i] = io_panArray[i_nSize - (1+i)];
    io_panArray[i_nSize - (1+i)] = nTmp;
  }
}


//==============================
// wrapper for metis
//==============================
void  REORD_NDOrdering::DoMetisND(
   REORD_QGraph *i_psGraph,
  int *o_panPerm,
  int *o_panInvPerm)
{

  //Check input.
  assert(i_psGraph && o_panPerm && o_panInvPerm);

  int nNumV = i_psGraph->GetNumV();
  int *__restrict panXAdj = i_psGraph->GetXAdj();
  assert(panXAdj);
  int *__restrict panAdj = i_psGraph->GetAdj();
  assert(panAdj);

  int *vwgt = NULL;
  int options[40];
  for(int i = 0; i < 40; i++)
    options[i] = -1;

  int nStatus = 0;

  nStatus = METIS_NodeND(
              &nNumV,
              panXAdj,
              panAdj, 
              vwgt, 
              options, 
              o_panPerm,
              o_panInvPerm);
 
  if (nStatus != 1) 
    printf("\n***Metis returned with an error.\n");
}

//+++++++++++++++++++++++++++++++++++++++++++
// Class Util
//
// Mehdi Paak
//+++++++++++++++++++++++++++++++++++++++++++
double Util::TimeDiff(double i_dTime)
{

  timeb TimeNow;
   REORD_me(&TimeNow);
  double Dt = (double)(TimeNow.time) + (double)(TimeNow.millitm)/1000.-i_dTime;
  return Dt;
}

void Util::ReadSparseMatCRS(
	const char *i_file, 
	int **o_ptr,  // IA       
	int **o_indx,  // JA
	double **o_val,  // A
	int *o_n, // Matrix nxn
	int *o_nnz, // Number of nonzeros
	const unsigned i_nBuffSize) 
{

	FILE * pFile;
  fopen_s(&pFile,i_file,"rb");

	if(!pFile)
	{
		fprintf(stderr, "Error opening file\n");
		return;
	}

	unsigned  nNumEq, nNNZ, i;
	
	// 1. Read number of equations.
	if(fread(&nNumEq, sizeof(nNumEq), 1, pFile) != 1)
	{
		fprintf(stderr, "Error reading number of equations \n");
		fclose(pFile);
		return;	// TODO: return int check in main for exit.
	}


	// Read data sparse matrix in crs format pointer (IA), index (JA), value (A); 0-based indexing 
	// 1. Poiter array (IA)

	//int *ptr = (int*) malloc((nNumEq + 1) * sizeof(*ptr));
  int *ptr = new int[nNumEq + 1];
	
	unsigned nPtrStps = (nNumEq+1)/i_nBuffSize;  
	unsigned nPtrRm =   (nNumEq+1)%i_nBuffSize;
	
    
    i=0;
	for(i = 0; i < nPtrStps * i_nBuffSize; i += i_nBuffSize)
	{
      if(fread(ptr+i, sizeof(*ptr), i_nBuffSize, pFile) != i_nBuffSize)
	  {
		fprintf(stderr, "Error while reading pointer array \n");
		fclose(pFile);
		return;	
	  } 
	}
  
	if(nPtrRm > 0)
	  if(fread(ptr+i, sizeof(*ptr), nPtrRm, pFile) != nPtrRm)
	  {
		fprintf(stderr, "Error while reading pointer array \n");
		fclose(pFile);
		return;	
	  } 		

	nNNZ = ptr[nNumEq];

	unsigned nIndStps = nNNZ/i_nBuffSize;
	unsigned nIndRm = nNNZ%i_nBuffSize;	
	
	// 2. Index array (JA)

	//int *indx = (int*) malloc(nNNZ * sizeof(*indx));
  int *indx = new int[nNNZ];
  
  i=0;
	for(i = 0; i < nIndStps * i_nBuffSize; i += i_nBuffSize)
	{
	  if(fread(indx+i, sizeof(*indx), i_nBuffSize, pFile) != i_nBuffSize)
	  {
		fprintf(stderr, "Error while reading index array \n");
		fclose(pFile);
		return;	
	  } 
	}

	if(nIndRm > 0)
	  if(fread(indx+i, sizeof(*indx), nIndRm, pFile) != nIndRm)
	  {
		fprintf(stderr, "Error while reading pointer array \n");
		fclose(pFile);
		return;	
	  } 
	  
	// 3. Value array (A)
  i=0;
	//double *val = (double*) malloc(nNNZ * sizeof(val));
  double *val = new double[nNNZ];

	for(i = 0; i < nIndStps * i_nBuffSize; i += i_nBuffSize)
	{
	  if(fread(val+i, sizeof(*val), i_nBuffSize, pFile) != i_nBuffSize)
	  {
	  	fprintf(stderr, "Error while reading value array \n");
		fclose(pFile);
		return;	
	  } 
    }
	
	if(nIndRm > 0)
	  if(fread(val+i, sizeof(*val), nIndRm, pFile) != nIndRm)
	  {
		fprintf(stderr, "Error while reading pointer array \n");
		fclose(pFile);
		return;	
	  } 	
	
	
	fclose(pFile);

	*o_n = nNumEq;
	*o_nnz = nNNZ;
  *o_ptr =ptr;
  *o_indx = indx;
  *o_val = val;
}

void Util::ReadRhs(
	 const char *i_file, // File name
   const int i_n,   // Size
	 double **o_RhsVec,
   const unsigned i_nBuffSize)
{
	FILE * pFile;
  fopen_s(&pFile,i_file,"rb");
	if(!pFile)
	{
		fprintf(stderr, "Error opening file\n");
		return;
	}

  // Chunk sizes.
	unsigned nStps = (i_n)/i_nBuffSize;  
	unsigned nRm =   (i_n)%i_nBuffSize;



	// Read 1-D array of size i_n; 0-based indexing 
	//double * RhsVec = (double*) malloc(i_n * sizeof(*RhsVec));
  double * RhsVec = new double[i_n];

	unsigned i;
  for(i = 0; i < nStps * i_nBuffSize; i += i_nBuffSize)
  	if(fread(RhsVec + i, sizeof(*RhsVec), i_nBuffSize, pFile) != i_nBuffSize )
  	{
	  	fprintf(stderr, "Error while reading from array \n");
	  	fclose(pFile);
	  	return;	
	  }

  if(nRm > 0)
  	if(fread(RhsVec + i, sizeof(*RhsVec), nRm, pFile) != nRm )
  	{
	  	fprintf(stderr, "Error while reading from array \n");
	  	fclose(pFile);
	  	return;	
	  }

	fclose(pFile); 

  *o_RhsVec = RhsVec;
}

void Util::PrintVecToFile(
	const int i_n, 
	const int *i_panVec, 
	const char *i_pcFileName)
{
  // Check input.
	assert(i_n > 0 && i_panVec && i_pcFileName);

  // Open file.
	FILE *pfFile;
	if (fopen_s(&pfFile, i_pcFileName, "w"))
		fprintf(stdout, "\nError creating file %s\n", i_pcFileName);

	for (int i = 0; i < i_n; i++)
	{
		fprintf(pfFile, "%d\n", i_panVec[i]);
	}
  // Close file.
	fclose(pfFile);
}