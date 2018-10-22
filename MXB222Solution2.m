%%Initialise the coefficient matrix using 5 point central differnence method. 
InitialiseVariables;

%%Check if coefficient matrix is SPD or not
isSPD(A);

%Generate the data structures
packed           = GenerateDataStructures("packed", A);
[banded, RCM]    = GenerateDataStructures("band", A);
[sprse , AMD]    = GenerateDataStructures("sparse", A);
[CSR,csr_row,csr_col]     = GenerateDataStructures("CSR", A);

%Generate cholesky factoristion of the data structures above. 
full_chol                           = chlsky("full", A);
[packed_chol,packed_chol_linear]    = chlsky("packed", packed);
banded_chol                         = chlsky("band", banded);
sprse_chol                         = chlsky("sparse", sprse);
CSR_chol                            = chlsky("CSR", CSR);

%Uncomment code below to generate all plots
plots;
