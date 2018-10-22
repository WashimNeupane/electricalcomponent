function [A, A_packed] = chlsky(type, A)

%%INPUT PARAMETERS
%matrix =>    Pass on a matrix in either packed, banded, sparse or CSR form. 
%type == F :  the cholesky factor for full   2d array. 
%type == P :  the cholesky factor for packed 1d array. 
%type == B :  the cholesky factor for band   2d array. 
%type == S :  the cholesky factor for sparse 2d array. 
%type == C :  the cholesky factor for CSR    2d array. 

if(type == "full")
    n = size(A);
    for j = 1:n 
      for i = 1:j-1        
        for k = 1:i-1 
            A(i,j) = A(i,j) - A(k,i) * A(k,j); 
         end              
        A(i,j) = A(i,j) / A(i,i);
        A(j,j) = A(j,j) - A(i,j)^2; 
       end
    A(j,j) = sqrt(A(j,j));
    end
 
elseif(type == "packed")   
    m = 33;
    A_packed = A;
    for j = 1:m             
      for i = 1:j-1         
         for k = 1:i-1                
             A_packed(i+(j*(j-1)/2)) = A_packed(i+(j*(j-1)/2)) - A_packed(k+(i*(i-1)/2)) * A_packed(k+(j*(j-1)/2)); 
         end              
         A_packed(i+(j*(j-1)/2)) = A_packed(i+(j*(j-1)/2)) / A_packed(i+(i*(i-1)/2));
         A_packed(j+(j*(j-1)/2)) = A_packed(j+(j*(j-1)/2)) - A_packed(i+(j*(j-1)/2))^2; 
       end
    A_packed(j+(j*(j-1)/2)) = sqrt(A_packed(j+(j*(j-1)/2)));
    end
    
 %%unpacks the factorisation above into a matrix
 
 A = zeros(33);
 for j= 1:m
    for i = 1:j       
         A(i,j)= A_packed(i+(j*(j-1)/2)); 
    end
 end
 
  elseif(type == "band")
[m,p] = size(A);  p = p - 1;
for k = 1:m
    last = min(k+p,m) - k + 1;
    for j = 2:last
        i = k + j - 1;
        A(i,1:last-j+1) = A(i,1:last-j+1) - ((A(k,j))/A(k,1))*A(k,j:last); 
    end
    A(k,:) = real(A(k,:)/sqrt(A(k,1)));
end
A(end-(p-1):end,end) = 0;
    
elseif(type == "sparse")
    %%input code here
    
elseif(type == "CSR")
    %%input code here
end
end
