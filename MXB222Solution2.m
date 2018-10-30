%%Initialise the coefficient matrix using 5 point central differnence method. 
[A,b,T] = InitialiseVariables;

%%Check if coefficient matrix is SPD or not
isSPD(A);

%Take upper part of A
% A = triu(A);

%Generate the data structures
full                             = triu(A);
packed                           = GenerateDataStructures("packed", A);
[banded, banded_RCM]             = GenerateDataStructures("band", A);
[sprse , sparse_AMD]             = GenerateDataStructures("sparse", A);
[CSR,CSR_RRow,CSR_RCol]          = GenerateDataStructures("CSR", A);

%%Generate cholesky factoristion of the data structures above. 
%%Use forward and backward substitution to generation numerical solution 
full_chol                        = chlsky("full", full);
z                                = substitution("forward", full_chol', b);
X_full                           = substitution("backward", full_chol, z); 

%packed_chol is full matrix, packed_chol_linear is packed matrix
[packed_chol,packed_chol_linear] = chlsky("packed", packed);
z                                = substitution("forward", packed_chol', b);
X_packed                         = substitution("backward", packed_chol, z); 

%Our own implementation of cholesky for banded matrices.
banded_chol                      = chlsky("band",  banded);

%%Using chol to generate an unpacked version. Comparing the above matrix
%%with the matrix generated below, we can say that they are identical. 
banded_chol_unpacked             = chol(banded_RCM);                    
z                                = substitution("forward", banded_chol_unpacked', b);
X_banded                         = substitution("backward", banded_chol_unpacked, z); 

sprse_chol                       = chlsky("sparse", triu(sprse));
z                                = substitution("forward", sprse_chol', b);
X_sprse                          = substitution("backward", sprse_chol, z); 

i = 11; %runtime
j = 5; %iterations
%%Generate the solution to CSR using iterative methods
[CSR_Jacobi, iteration(j,1),runtime(i,1)]                  = IterativeMethodCSR("Jacobi", CSR,CSR_RRow,CSR_RCol,b);
[CSR_GaussSeidel, iteration(j,2),runtime(i,2)]             = IterativeMethodCSR("GaussSiedel", CSR,CSR_RRow,CSR_RCol,b);
[CSR_ConjugateGradient, iteration(j,3),runtime(i,3)]       = IterativeMethodCSR("ConjugateGradient", CSR,CSR_RRow,CSR_RCol,b);

%Generate SOR with varying omega values
[CSR_SOR1,iteration(j,4),runtime(i,4)]          = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.1);
[CSR_SOR2,iteration(j,5),runtime(i,5)]          = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.5);
[CSR_SOR3,iteration(j,6),runtime(i,6)]          = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.94);
[CSR_SOR4,iteration(j,7),runtime(i,7)]          = IterativeMethodCSR("SOR", CSR,CSR_RRow,CSR_RCol,b,0.77);

%Solutions using matlab's library. This output is compared against other to
%check for errors
Sol = A\b; 

%Uncomment code below to generate all plots
plots;
efficiency;
