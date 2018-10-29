%%Initialise the coefficient matrix using 5 point central differnence method. 
[A,b,T] = InitialiseVariables;

%%Check if coefficient matrix is SPD or not
isSPD(A);

%Take upper part of A
% A = triu(A);

%Generate the data structures
full                         = triu(A);
packed                       = GenerateDataStructures("packed", A);
[banded, banded_RCM]         = GenerateDataStructures("band", A);
[sprse , sparse_AMD]         = GenerateDataStructures("sparse", A);
[CSR,CSR_RRow,CSR_RCol]      = GenerateDataStructures("CSR", A);

%%Generate cholesky factoristion of the data structures above. 
full_chol                           = chlsky("full", full);

%Note: packed_chol is full matrix, packed_chol_linear is packed matrix
[packed_chol, packed_chol_linear]   = chlsky("packed", packed); 
banded_chol                         = chlsky("band", banded);
sprse_chol                          = chlsky("sparse", sprse);

%Generate the solution to CSR using iterative methods
[CSR_Jacobi, CSR_Jacobi_numIterations]                       = IterativeMethodCSR("Jacobi", CSR,CSR_RRow,CSR_RCol,b);
[CSR_GaussSeidel, CSR_GaussSeidel_numIterations]             = IterativeMethodCSR("GaussSiedel", CSR,CSR_RRow,CSR_RCol,b);
[CSR_ConjugateGradient, CSR_ConjugateGradient_numIterations] = IterativeMethodCSR("ConjugateGradient", CSR,CSR_RRow,CSR_RCol,b);

%Generate SOR with varying omega values
[CSR_SOR1,SOR1_numIterations]                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.01);
[CSR_SOR2,SOR2_numIterations]                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.5);
[CSR_SOR3,SOR3_numIterations]                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.94);
[CSR_SOR4,SOR4_numIterations]                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.77);

%Solutions using direct methods
Sol_full = full\b; 

%Uncomment code below to generate all plots
plots;
