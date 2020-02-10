//-----------------------------------------------
// Class for graph manipulations and tools for sparse matrix 
// Reordering. The Reordering is performed to minimize the
// number of nonzeros in the factorized matrix.
// Algorithms are based on the reordering algorithms 
// proposed by Alan George and  JWH Liu.
// See https://epubs.siam.org/doi/abs/10.1137/1031001
//
// Minimum Degree (MD), Multiple minimum degree (MMD), Nested Disection (ND)
//
// By: Mehdi Paak
// TODO: major refactoring needed
//----------------------------------------------

#pragma once
#define MinDeg_EXPORTS

#ifdef MinDeg_EXPORTS
#define MinDegLIBRARY_API __declspec(dllexport) 
#else
#define MinDegLIBRARY_API __declspec(dllimport) 
#endif

#define SDS_INT __int32 //__int64 : long long

enum Indexing 
{
  C_Style,
  F_Style
};

enum ReorderingMethod
{
  nOrig,
  nMD,
  nMMD,
  nND,
  nMetisND
};

enum MatForm
{
  nUpper,
  nLower,
  nFull
};

//METIS_NodeND(int nvtxs, int *xadj, int *adjncy, int *vwgt, int *options, int *perm,int *iperm);
// Metis implementation of Nested Disection
extern "C" int METIS_NodeND(int *nvtxs, int *xadj, int *adjncy, int *vwgt, int *options, int *perm,int *iperm);

class MinDegLIBRARY_API  REORD_QGraph
{

public:

  // Constructor.
   REORD_QGraph():m_panXAdj(0), m_panAdj(0), m_nNumV(0),m_nNumE(0), m_nStyle(C_Style){}; 
   REORD_QGraph(const int i_nNumV, const int i_nNumE, Indexing i_Style = C_Style);
  
  // Copy constructor.
   REORD_QGraph(const  REORD_QGraph &i_oGraph);

  // Destructor.
  ~ REORD_QGraph();

  void SetDataSpace(
    const int i_nNumV, 
    const int i_nNumE);

  //===== Get infos ======.

  // Get vertex number.
  int GetNumV() const {return m_nNumV;}

  // Get Edge number.
  int GetNumE() const {return m_nNumE;}

  // Get XAdj.
  const int *GetXAdj() const {return m_panXAdj;}
  int *GetXAdj() {return m_panXAdj;}
  
  // Get Adj.
  const int *GetAdj() const {return m_panAdj;}
  int *GetAdj() {return m_panAdj;}
  
  // Is empty
  const bool IsEmpty() const {return (m_nNumV == 0);}
  
  // Get Style
  const Indexing GetStyle() const {return m_nStyle;}
  Indexing GetStyle() {return m_nStyle;}
  
  // Print Graph: XAdj \n Adj
  void PrintGraph() const;
  void PrintToFile(const char *i_pcFileName);
  //===== Set data =======.
  void SetNumVTo(const int i_nNumV) {m_nNumV = i_nNumV;}
  void SetNumETo(const int i_nNumE) {m_nNumV = i_nNumE;}
  void SetStyleTo(const Indexing i_Style){m_nStyle = i_Style;}

  // Convert to Fortran style 1-based indexing
  void ConvertToFort();
 
  // Convert from Fort to C style 0-based indexing
  void ConvertToC();

  // Read graph from file.
static void ReadGrphFromFile(
   REORD_QGraph *i_psQGraph,
  const char *i_pacFileName);

static void ChangeToQGrph(
  const int &i_nNumV,
  const int &i_nNumE,
  const char *i_pacInFileName,
  const char *i_pacOutFileName);

protected:

  int *m_panXAdj;
  int *m_panAdj;
  int m_nNumV;       // Number of vertices, i.e., number of equations
  int m_nNumE;       // Number of edges, i.e., number of nonzeros - number of diagonals (NNZ - nNumV) excluding self-edges.
  Indexing m_nStyle;
};

//=========================================
class SparseMat:public  REORD_QGraph
{

public:
  // Constructor.
  SparseMat():m_padVal(0){}; //  REORD_QGraph() is called automatically first.
  SparseMat(const int i_nNumV, const int i_nNumE, Indexing i_Style = C_Style);
  
  // Copy constructor.
  SparseMat(const SparseMat &i_psSparseMat);

  // Destructor.
  ~SparseMat();

  // Allocate space.
  void SetDataSpaceSparseMat( 
    const int i_nNumE);

  // Get value.
  const double *GetVal() const {return m_padVal;}
  double *GetVal() {return m_padVal;}

  // Set Values.
  void SetSpMatInfo(const int &i_nNumEq, const int &i_nNNZ, int *i_panPtr, int *i_panIndx, double *i_padVal)
  {
    m_nNumV = i_nNumEq;
    m_nNumE = i_nNNZ;
    m_panXAdj = i_panPtr;
    m_panAdj = i_panIndx;
    m_padVal = i_padVal;
  }

  // Form graph adjaceny from matrix upper triangular part ordered.
  static void FormGraphAdj(const SparseMat *i_psU,  REORD_QGraph *o_psQGraph, int &o_nMaxRowNNZ);

  // Get Diagonal.
  static void GetDiagonal(const SparseMat *i_psU, const MatForm i_nMatForm, double *o_padDiag);

  // Insert A into L (NNZL>>NNZA) based on permutation vec, Diagonals are stored separately.
  static void InsertAToL(
    const SparseMat *i_psU, 
    const MatForm i_nMatForm, 
    const int *i_panPerm, 
    const int *i_panInvPerm, 
    double *o_padLnz,
    double *o_padDiag);


private:
  double *m_padVal;
};

//==============================================================
// MinDeg data structure
// F_Style
//==============================================================
struct MinDegLIBRARY_API  REORD_MinDegData_s
{
 // Constructor.
   REORD_MinDegData_s(): m_panDeg(0),m_panMarker(0), m_panQLink(0), m_panQSize(0), m_nNumV(0),m_nNumE(0),
    m_panRchSet(0), m_panNbrhdSet(0), m_nRchSize(0), m_nNbrhdSize(0), m_nRchTotNodes(0),m_nNumElimSupNodes(0),m_nStyle(C_Style){}; 

   REORD_MinDegData_s(const int i_nNumV); // create C_Style data arrays. 
 // Destructor.
  ~ REORD_MinDegData_s();

 //====== Get info ==========

  // Get vertex number.
  int GetNumV() const {return m_nNumV;}
 
  // Get Deg.
  const int *GetDeg() const {return m_panDeg;}
  int *GetDeg() {return m_panDeg;}
  
  // Get Marker.
  const int *GetMarker() const {return m_panMarker;}
  int *GetMarker() {return m_panMarker;}
  
  // Get QLink.
  const int *GetQLink() const {return m_panQLink;}
  int *GetQLink() {return m_panQLink;}

  // Get QSize.
  const int *GetQSize() const {return m_panQSize;}
  int *GetQSize() {return m_panQSize;}


  // Is empty.
  const bool IsEmpty() const {return (m_nNumV == 0);}
  
  // Get Rch set size.
  const int GetRchSize() const {return m_nRchSize;}

  // Get Neighborhood set size.
  const int GetNbrhdSize() const {return m_nNbrhdSize;}

  // Get Rch set total nodes. 
  const int GetRchTotNodesNum() const {return m_nRchTotNodes;}

  // Get number of eliminated supernodes.
  const int GetNumElimSupNodes() const {return m_nNumElimSupNodes;}
  
  // Get Reach Set.
  const int *GetRchSet() const {return m_panRchSet;}
  int *GetRchSet() {return m_panRchSet;}

  // Get Nbrhd Set.
  const int *GetNbrhdSet() const {return  m_panNbrhdSet;}
  int *GetNbrhdSet() {return  m_panNbrhdSet;}


  // Get Style
  const Indexing GetStyle() const {return m_nStyle;}
  Indexing GetStyle() {return m_nStyle;}

  static void PrintVec(
    const int &i_nLen, 
    const int *i_panVec, 
    const Indexing i_Style = C_Style);

  //====== Set methods ==========

  void SetStyleTo(const Indexing i_Style){m_nStyle = i_Style;}

  void SetRchSizeTo(int &i_nRchSize){m_nRchSize = i_nRchSize;}
  void SetNbrhdSizeTo(int &i_nNbrhdSize){m_nNbrhdSize = i_nNbrhdSize;}
  
  void SetRchTotNodesNumTo(int &i_nRchTotNodesNum) {m_nRchTotNodes = i_nRchTotNodesNum;}
  void SetNumElimSupNodesTo(int &i_nNumElimSupNodes) {m_nNumElimSupNodes = i_nNumElimSupNodes;}

  void ConvertToFort();
  void ConvertToC();

private:

  int *m_panDeg;     // Array containing the degree of each vertex.
  int *m_panMarker;  // A work array heavily used in MinDeg algs: -1 means vertices joined in a super node etc.
  int *m_panQLink;   // Linked list for indistinguishable supernodes.
  int *m_panQSize;   // Integer array indicating supernode size associated with the supernode head.

  int *m_panRchSet;
  int *m_panNbrhdSet;
  int  m_nRchSize;
  int  m_nNbrhdSize;
  int  m_nRchTotNodes; // i.e. nRchSize + additional nodes in Xsupernodes (nodes in Rch set which are indistinguishable.). 
  int  m_nNumElimSupNodes; // Number of Eliminated supernodes adjacent to some node in RchSet.

  int m_nNumV;       // Number of vertices, i.e., number of equations
  int m_nNumE;       // Number of edges, i.e., number of nonzeros - number of diagonals (NNZ - nNumV) excluding self-edges.

  Indexing m_nStyle;
};

//==============================================================
// Static class for SDS operations on the linear system graph.
// 
// Mehdi P
//==============================================================
class MinDegLIBRARY_API  REORD_SDSMinDeg
{
  
public:
  
//===========================
// Data initialization.
//===========================
  static void InitDataForMD(
     REORD_QGraph *io_psQGraph, 
     REORD_MinDegData_s *o_psMinDegData,
    int * __restrict o_panPerm,
    int * __restrict o_panInvPerm,
    int &o_nFistMinDeg);

//===========================
// Find Min Deg permutation.
//===========================
  static void FindMinDegPerm(
     REORD_QGraph *io_psQGraph,
    int * __restrict o_panPerm,
    int * __restrict o_panInvPerm,
    SDS_INT & o_nNumSubs);

//===========================
// Find Min Deg permutation.
//===========================
  static void FindMinDegPerm_2(
     REORD_QGraph *io_psQGraph,
    int * __restrict o_panPerm);

//==========================
// Static function to find the Reachable set from the node Root
// through the eliminated set S.
// Reach(y,S) = Adj(Nbrhd(y,S) U {y})
// works with Fortran style 1-indexing
//==========================
  static void FindReachSet(
    const int i_nVertex,
     REORD_QGraph *io_psQGraph,
     REORD_MinDegData_s *o_psMinDegData);

//==========================
// Static function to find the Reachable set from the node Root
// through the eliminated set S.
// Reach(y,S) = Adj(Nbrhd(y,S) U {y})
// works with Fortran style 1-indexing
//==========================
  static void FindReachSet(
    const int i_nVertex,
     REORD_QGraph *io_psQGraph,
    int &o_nRchZise,
    int *__restrict o_panRchSet,
    int &o_nNbrhdSize,
    int *__restrict o_panNbrhdSet,
     REORD_MinDegData_s *o_psMinDegData);

//=============================
// static function for the quotient graph transformation after 
// elimination of a node or xnode
// First Reachset is added to the representative of new supernode,
// Then this representative is added to the adj list of Reachset nodes.
// works with Fortran style 1-indexing.
//=============================
  static void TransformQuotientGraph(
    const int i_nVertex,
     REORD_QGraph *__restrict io_psQGraph,
    const  REORD_MinDegData_s * __restrict i_psMinDegData);

//============================
// static function for updating degrees and merging indistinguishable supernodes.
//
// Works with F_Style.
//============================
  static void UpdateDeg(
     REORD_QGraph * __restrict io_psQGraph,
     REORD_MinDegData_s * i_psMinDegData);

//============================
// static function for finding and merging indistinguishable supernodes (XNodes).
//
// Works with F_Style.
//============================
  static void MergeXNodes(
     REORD_QGraph * __restrict io_psQGraph,
     REORD_MinDegData_s * __restrict i_psMinDegData);

};

//==============================================================
// MinDeg data structure
// F_Style
//==============================================================
struct MinDegLIBRARY_API  REORD_MMD_s
{
 // Constructor.
   REORD_MMD_s(): m_panHead(0),m_panFwd(0), m_panBwd(0), m_panQsize(0), m_panList(0), 
    m_panMarker(0), m_nNumV(0),m_nStyle(C_Style){}; 

   REORD_MMD_s(const int i_nNumV); // create C_Style data arrays. 
 // Destructor.
  ~ REORD_MMD_s();

 //====== Get info ==========

  // Get vertex number.
  int GetNumV() const {return m_nNumV;}
 
  // Get Degree Head list.
  const int *GetHead() const {return m_panHead;}
  int *GetHead() {return m_panHead;}
  
  // Get Forward list.
  const int *GetFwd() const {return m_panFwd;}
  int *GetFwd() {return m_panFwd;}
  
  // Get Backward list.
  const int *GetBwd() const {return m_panBwd;}
  int *GetBwd() {return m_panBwd;}

  // Get QSize.
  const int *GetQsize() const {return m_panQsize;}
  int *GetQsize() {return m_panQsize;}

  // Get List.
  const int *GetList() const {return m_panList;}
  int *GetList() {return m_panList;}

  // Get Marker.
  const int *GetMarker() const {return m_panMarker;}
  int *GetMarker() {return m_panMarker;}

  // Is empty.
  const bool IsEmpty() const {return (m_nNumV == 0);}
  
  // Get Style
  const Indexing GetStyle() const {return m_nStyle;}
  Indexing GetStyle() {return m_nStyle;}

  //====== Set methods ==========

  void SetStyleTo(const Indexing i_Style){m_nStyle = i_Style;}
  void SetFwdBckwdPointers (int *i_panInvPerm, int *i_panPerm);
  void ConvertToFort();
  void ConvertToC();

private:

  int *m_panHead;     // Deg Head for doubly-linked list Head Forward Bakward.
  int *m_panFwd;      // Deg Forward list.
  int *m_panBwd;      // Deg Backward list.
  int *m_panQsize;    // List for XNodes (indistinguishable).
  int *m_panList;     // Temporary linked list.
  int *m_panMarker;   // Marker array.

  int m_nNumV;       // Number of vertices, i.e., number of equations.

  Indexing m_nStyle;
};

//==============================================================
// Static class for MMD.
// F_Style
// Mehdi P
//==============================================================
class MinDegLIBRARY_API  REORD_SDSMMD
{

public:
//===========================
// Find MMD permutation.
//===========================
  static void FindMMDPerm(
    const int &i_nDelta,
     REORD_QGraph *io_psQGraph,
    int *  o_panPerm,
    int *  o_panInvPerm,
    SDS_INT &o_nNumSubs);


//===========================
// Eliminate MinDeg node and transform quotient graph.
//===========================
  static void MMDQGraphTransform(
    const int &i_nMDegVtx,
     REORD_QGraph *io_psQGraph,
     REORD_MMD_s *io_psMMDData,
    const int i_nMXINT,
    const int i_nTag);

//===========================
// Update degree of updatable nodes.
//===========================
  static void MMDUpdateDeg(
    const  REORD_QGraph *i_psQGraph,
     REORD_MMD_s *io_psMMDData,
    const int i_nEHead,
    const int i_nDelta,
    int &io_nMinDeg,
    const int i_nMXINT,
    int &io_nTag);

//===========================
// Data initialization.
//===========================
  static void InitDataMMD(
    const  REORD_QGraph *i_psQGraph, 
     REORD_MMD_s *o_psMMDData);


//===========================
// Clculate the perm and inv perm 
// numbering from info already contained in them.
//===========================
  static void MMDNumber(
    const int i_nNumV,
    const  REORD_MMD_s *i_psMMDData,
    int * io_panInvPerm,
    int * io_panPerm);





};

//==============================================================
// Static class for ND.
// Mehdi P
//==============================================================
class MinDegLIBRARY_API  REORD_NDOrdering
{
public:

//===========================
// Find Level structure rooted at Root. (BFS)
//===========================
  static void FindLevelStruct(
    const int &i_nRoot,
    const int *__restrict i_panXAdj,
    const int *__restrict i_panAdj,
    int *__restrict io_panMask,
    int &o_nNumLvl,
    int *__restrict o_panXLS,
    int *__restrict o_panLS);

//===========================
// Find Pseudo-peripheral node.
//===========================
  static void PsedoPeripheralNodeLS(
    int &io_nRoot,
    const int *__restrict i_panXAdj,
    const int *__restrict i_panAdj,
    int *__restrict i_panMask,
    int &o_nNumLvl,
    int *__restrict o_panXLS,
    int *__restrict o_panLS);

//===========================
// Nested Diesection permutation.
//===========================
  static void DoNDOrdering(
     REORD_QGraph *i_psGraph,
    int *o_panPerm,
    int *o_panInvPerm);

//===========================
// Reverse integer array.
//===========================
  static void ReverseArray(
    const int i_nSize,
    int *io_panArray);

//===========================
// Find Separator and number.
//===========================
  static void FindSeparator(
    int &i_nRoot,
    int *__restrict i_panXAdj,
    const int *__restrict i_panAdj,
    int *__restrict io_panMask,
    int &o_nSepSize,
    int *__restrict o_panSepVec,
    int *__restrict panXLS,
    int *__restrict panLS);

//==============================
// wrapper for metis
//==============================
  static void DoMetisND(
     REORD_QGraph *i_psGraph,
    int *o_panPerm,
    int *o_panInvPerm);


};

//==============================================================
//
//
//
//==============================================================
class MinDegLIBRARY_API  REORD_Solve
{

public:

  static void DoSymbFact(
    const  REORD_QGraph *i_psGraph,
    const int *i_panPerm,
    const int *i_panInvPerm,
    int *o_panXLnz,
    int *o_panXNzsub,
    int *o_panNzSub,
    int &o_nNumNNZ,
    int &o_nNumSubs,
    int &o_eError);

  static void DoSymbFact2(
    const  REORD_QGraph *i_psGraph,
    const int *i_panPerm,
    const int *i_panInvPerm,
    SDS_INT *o_panXLnz,
    SDS_INT *o_panXNzsub,
    SDS_INT **o_panNzSub,
    SDS_INT &o_nNumNNZ,
    SDS_INT &o_nNumSubs,
    int &o_eError);

};

//=======================================
class Util
{

public:
  
  static double TimeDiff(double i_dTime);
  
  static void ReadSparseMatCRS(
	  const char *i_file, 
	  int **o_ptr,  // IA       
	  int **o_indx,  // JA
	  double **o_val,  // A
	  int *o_n, // Matrix nxn
	  int *o_nnz, // Number of nonzeros
	  const unsigned i_nBuffSize = 1024 * 1024);

  static  void ReadRhs(
	  const char *i_file, // File name
    const int i_n,   // Size
	  double **o_RhsVec,
    const unsigned i_nBuffSize = 1024);

  static void PrintVecToFile(const int i_n, const int *i_panVec, const char *i_pcFileName);
};

