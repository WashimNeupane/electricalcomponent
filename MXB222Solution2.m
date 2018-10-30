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
%%Use forward and backward substitution to generation numerical solution 
full_chol                        = chlsky("full", full);
z = substitution("forward", full_chol', b);
X_full = substitution("backward", full_chol, z); 

%packed_chol is full matrix, packed_chol_linear is packed matrix
[packed_chol,packed_chol_linear] = chlsky("packed", packed);
z = substitution("forward", packed_chol', b);
X_packed = substitution("backward", packed_chol, z); 

%banded_chol_unpacked is full matrix, banded_chol is banded
banded_chol_unpacked = chlsky("full",   triu(banded_RCM));
banded_chol                      = GenerateDataStructures("band",banded_chol_unpacked);
z = substitution("forward", banded_chol_unpacked', b);
X_banded = substitution("backward", banded_chol_unpacked, z); 

sprse_chol                       = chlsky("sparse", triu(sprse));
z = substitution("forward", sprse_chol', b);
X_sprse = substitution("backward", sprse_chol, z); 

%%Generate the solution to CSR using iterative methods
[CSR_Jacobi, iteration(1)]                  = IterativeMethodCSR("Jacobi", CSR,CSR_RRow,CSR_RCol,b);
[CSR_GaussSeidel, iteration(2)]             = IterativeMethodCSR("GaussSiedel", CSR,CSR_RRow,CSR_RCol,b);
[CSR_ConjugateGradient, iteration(3)]       = IterativeMethodCSR("ConjugateGradient", CSR,CSR_RRow,CSR_RCol,b);

%Generate SOR with varying omega values
[CSR_SOR1,iteration(4)]                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.1);
[CSR_SOR2,iteration(5)]                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.5);
[CSR_SOR3,iteration(6)]                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.94);
[CSR_SOR4,iteration(7)]                            = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.77);

%Solutions using direct methods
Sol_full = A\b; 

%Uncomment code below to generate all plots
plots;
