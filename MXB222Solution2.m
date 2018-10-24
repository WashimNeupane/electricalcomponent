%%Initialise the coefficient matrix using 5 point central differnence method. 
[A,b,T] = InitialiseVariables;

%%Check if coefficient matrix is SPD or not
isSPD(A);

%Generate the data structures
packed                       = GenerateDataStructures("packed", A);
[banded, banded_RCM]         = GenerateDataStructures("band", A);
[sprse , sparse_AMD]         = GenerateDataStructures("sparse", A);
[CSR,CSR_RRow,CSR_RCol]      = GenerateDataStructures("CSR", A);

%Generate cholesky factoristion of the data structures above. 
full_chol                           = chlsky("full", A);
[packed_chol,packed_chol_linear]    = chlsky("packed", packed);
banded_chol                         = chlsky("band", banded);
sprse_chol                          = chlsky("sparse", sprse);

%Generate the solution to CSR using iterative methods
CSR_Jacobi                          = IterativeMethodCSR("Jacobi", CSR,CSR_RRow,CSR_RCol,b);
CSR_GaussSeidel                     = IterativeMethodCSR("GaussSiedel", CSR,CSR_RRow,CSR_RCol,b);
CSR_ConjugateGradient               = IterativeMethodCSR("ConjugateGradient", CSR,CSR_RRow,CSR_RCol,b);

%Generate SOR with varying omega values
CSR_SOR1                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.01);
CSR_SOR2                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.5);
CSR_SOR3                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.94);
CSR_SOR4                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.77);

%Uncomment code below to generate all plots
plots;
