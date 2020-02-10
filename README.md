# Sparse-Matrix-Reordering
A comparison between several reordering algorithms for minimal non-zero fill-ins in the Cholesky factorization when solving Ax=b. 

The reordering algorithms implemented are:

- Minimum Degree (MD)
- multiple Minimum Degree (MMD)
- Nested Dissection (ND)
- Metis library Nested Dissection

Timing and the number of non-zeros in the factored matrix is reported.

# Usage: 
 MD.exe -filemat=<name> -filerhs=<name> -reordering=<0/1/2/3/4> 
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

 By: Mehdi Paak