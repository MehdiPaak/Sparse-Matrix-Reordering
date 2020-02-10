/*---------------------------------------------------------------------
 MD.cpp : Defines the entry point for the console application.
 Usage: MD.exe -filemat=<name> -filerhs=<name> -reordering=<0/1/2/3/4> 
 Example data files are provided in the sample data.
 
 -filemat: name of the binary file containing coeeficient matrix in 
 compressed row storage (CRS) format, the number of rows/equations is specified first
 then the pointer vector, column index vector and finally values
 
 -filerhs: file containing the right hand side to the equation A x= b (not necessary for reordering)
  
 
 -reordering: 
	0: no reordring
	1: Minimum Degree (MD)
	2: Multiple minimum degree (MMD)
	3: Nested Disectipon (ND)
	4: Nested Disection Metis library implementation
	
 - compare the number of nonzeros in the factorized matrix for comparing 
   the efficiency.
   
 Based on the reordering algorithms 
 proposed by Alan George and  JWH Liu.
 See https://epubs.siam.org/doi/abs/10.1137/1031001

 BY: Mehdi Paak
//--------------------------------------------------------------------*/

#include "stdafx.h"
#include "MinDeg.h"


int main(int argc, char* argv[])
{

  char* Mat_File_Name = 0;    
  char* Rhs_File_Name = 0; 
  int nReordering = 0;   // 0:MD, 1:MMD, 2:ND, 3:Metis
  int Reserved = 0;
  int k = 1;

  fprintf(stdout,"argc = %d\n",argc);

  if(argc > 1)
    while (k < argc) 
    {
      if ( !strncmp(argv[k],"-filemat",8) ) 
      {
        Mat_File_Name = &argv[k][9];
        k++;
      }  
      else if ( !strncmp(argv[k],"-filerhs",8 ) ) 
      {
  	    Rhs_File_Name = &argv[k][9];
        k++;
      }
      else if ( !strncmp(argv[k],"-reordering",11) )
      {
   	    nReordering = atoi(&argv[k][12]);
        k++;
      }  
      else if ( !strncmp(argv[k],"-spd", 4) )
      {
   	    Reserved = atoi(&argv[k][5]);
        k++;
      }  
      else 
      {
        fprintf(stderr, "invalide argument!\n");
        fprintf(stderr, "Usage: ./MD.exe -filemat=<name> -filerhs=<name> -reordering=<0/1/2/3/4> -spd=<0/1>\n");
        return -1;
      }
    }
  else
  {
    Mat_File_Name = "Mat_SB_Num.dat"; //"Mat_SB_Reinf_Roof.dat"; //"Mat_SB_Num.dat";
    Rhs_File_Name = "Rhs_SB_Num.dat"; //"Rhs_SB_Reinf_Roof.dat"; //"Rhs_SB_Num.dat";
    nReordering = 2;
  }



  printf("\n"
         "Mat files  =  %s\n"
         "Rhs File   =  %s\n"
         "Reordering = %d\n", Mat_File_Name, Rhs_File_Name, nReordering);
  

  // Matrix A, RHS and solution vectors.
  int nNumEq;
  int nNNZ;
  SparseMat sA;  // Ax = b.
  int *panPtr = NULL;
  int *panIndx = NULL; 
  double *padVal = NULL;
  double *padRhs = NULL;
  double *padSol = NULL;

  // Read Coefficient Matrix Sorted Upper Triangular including Diagonal.
  printf("\nRead Coeff. Matrix from %s ...\n", Mat_File_Name);
  double dT0 = Util::TimeDiff(0.0);

  Util::ReadSparseMatCRS(
    Mat_File_Name,
    &panPtr,
    &panIndx,
    &padVal,
    &nNumEq,
    &nNNZ);

  double dT1 = Util::TimeDiff(dT0);
  printf("Matrix read in %10.3f s\n", dT1);

  // Read Rhs.
  printf("\nRead Rhs from %s ...\n", Rhs_File_Name);
  dT0 = Util::TimeDiff(0.0);

  Util::ReadRhs(
    Rhs_File_Name,
    nNumEq,
    &padRhs);

  double dT2 = Util::TimeDiff(dT0);
  printf("Rhs read in %10.3f s\n", dT2);

  printf("\n"
         "NumEq = %d\n"
         "A NNZ = %d\n",nNumEq, nNNZ);

  sA.SetSpMatInfo(
    nNumEq,
    nNNZ,
    panPtr,
    panIndx,
    padVal);

  // Form graph adjacency list from U. 
  fprintf(stdout,"\nForming Adjacency list form U ...\n");
  int nMaxRowNNZ;
  dT0 = Util::TimeDiff(0.0);
  int nNumgraphE = 2 * (nNNZ - nNumEq);
   REORD_QGraph *psGraph = new  REORD_QGraph(nNumEq, nNumgraphE);
  SparseMat::FormGraphAdj(
    &sA, 
    psGraph, 
    nMaxRowNNZ);
  double dT3 =  Util::TimeDiff(dT0);
  fprintf(stdout,"Forming Adjacency done in %10.3f s\n",dT3);
  fprintf(stdout,"Row Max NNZ = %d\n\n", nMaxRowNNZ); 

  //++++++++ Reordering +++++++++++++++
  int nNumV = nNumEq;
  int *panPerm = new int [nNumV];
  int *panInvPerm = new int [nNumV];

  // Initialize.
  for(int i = 0; i < nNumV; i++)
  {
    panPerm[i] = i;
    panInvPerm[i] = i;
  }

  // Upperbound on the number of compressed subscripts.
  SDS_INT nNumSubs = 0; 
  
  dT0 = Util::TimeDiff(0.0);
   REORD_QGraph *psGraph_1 = NULL;
  switch(nReordering)
  {
    case nOrig:
      psGraph_1 = psGraph;
      nNumSubs = 5 * nNumV;
      fprintf(stdout,"\nOriginal Order, %d ...\n",nReordering);
      break;
  
    case nMD:
      // Destroys original graph copy for Symb. Fact.
      fprintf(stdout,"\nPerforming Reordering MD, %d ...\n",nReordering);
	    psGraph_1 =  new  REORD_QGraph (*psGraph);
       REORD_SDSMinDeg::FindMinDegPerm(
        psGraph,
        panPerm,
        panInvPerm,
        nNumSubs);
	    delete psGraph; 
      break;

    case nMMD:
	    // Destroys original graph copy for Symb. Fact.
      fprintf(stdout,"\nPerforming Reordering MMD, %d ...\n",nReordering);
	    psGraph_1 = new  REORD_QGraph(*psGraph);
       REORD_SDSMMD::FindMMDPerm(
        0,
        psGraph,     
        panPerm, 
        panInvPerm,
        nNumSubs); 
	    delete psGraph;
      break;

    case nND:
      // Innocuous.
      fprintf(stdout,"\nPerforming Reordering ND, %d ...\n",nReordering);
	    psGraph_1 = psGraph;
       REORD_NDOrdering::DoNDOrdering(
        psGraph,
        panPerm,
        panInvPerm);
      nNumSubs = 5 * nNumV; // In Symbfact Reallocates XNzSubs if needed.
      break;

    case nMetisND:
      fprintf(stdout,"\nPerforming Reordering MetisND, %d ...\n",nReordering);
      psGraph_1 = psGraph;
       REORD_NDOrdering::DoMetisND(
        psGraph,
        panPerm,
        panInvPerm);
      nNumSubs = 8 * (nNumV+nNumgraphE);
      break;

    default:
    fprintf(stderr,"Reordeing method %d is undefined (Orig:0, MD:1, MMD:2, ND:3, Metis: 4)", nReordering);
  }

  double dT4 = Util::TimeDiff(dT0);
  fprintf(stdout,"Reordering done in %10.3f s\n",dT4);
  fprintf(stdout,"NumSubs estimate %d \n",nNumSubs);

  //+++++++++ Symbolic Factorization +++++++++++++++++++++++++++++

  fprintf(stdout, "\nPerforming Symbolic Factorization ...\n");
  dT0 = Util::TimeDiff(0.0);
  SDS_INT nNNZFact = 0;
  int nError = 0;
  SDS_INT *panXLnz   = (SDS_INT*) malloc(sizeof(*panXLnz) * (nNumV + 1));
  SDS_INT *panXNzSub = (SDS_INT*)malloc(sizeof(*panXNzSub) * (nNumV + 1));
  SDS_INT *panNzSub  = (SDS_INT*)malloc(sizeof(*panNzSub) * nNumSubs);

   REORD_Solve::DoSymbFact2(
	  psGraph_1,
	  panPerm,
	  panInvPerm,
	  panXLnz,
	  panXNzSub,
	  &panNzSub,
	  nNNZFact,
	  nNumSubs, // io in case not enough memory.
	  nError);

  double dT5 = Util::TimeDiff(dT0);
  fprintf(stdout, "Symb. Fact. done in %10.3f s\n", dT5);
  fprintf(stdout, "NumSubs    =  %lld\n"
	                "Factor NNZ =  %lld\n\n", nNumSubs, nNNZFact);


  //+++++++++ Numerical Factorization +++++++++++++++++++++++++++++


 // Util::PrintVecToFile(nNumV, panPerm, "prm.txt");
 // sGraph.PrintToFile("G.txt");

//+++++++++++ Free Memory +++++++++++++++++++++

  delete [] panPerm;
  delete [] panInvPerm;
  delete [] padRhs;
  delete [] padSol;
  delete psGraph_1;

  free(panXLnz);
  free(panXNzSub);
  free(panNzSub);

  return 0;
}

