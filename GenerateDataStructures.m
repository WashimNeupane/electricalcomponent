function [A,B,C] = GenerateDataStructures(type, A)

if(type == "packed")
%Generating a packed matrix
[m,n] = size(A);
array = zeros(1, m+(n*(n-1)/2) );
for j= 1:m
    for i = 1:j       
          array(i+(j*(j-1)/2))= A(i,j);
    end
end
A = array;

elseif(type== "band")
%Generating a band matrix with RCM ordering
p = symrcm(A);
B = A(p,p);
A = triu(B);
A = spdiags(A);

elseif(type =="sparse")
%Generate a sparse matrix with AMD ordering
s = symamd(A);
B = A(s,s);
A = sparse(B);

elseif(type=="CSR")
%Generate a CSR ordered matrix
row= zeros(0);
col= zeros(0);
CSR= zeros(0);
for i=1:size(A)
    for j=1:size(A)
        if(A(i,j) ~= 0)
            row = [row ,i];
            col = [col,j];  
            CSR= [CSR, A(i,j)];
        end
    end
end

row_ord = ones(1);
for i=1:length(row)-1
    if(row(i+1)~= row(i))
        row_ord =[row_ord, i+1];     
    end
end
A = CSR;
B = row_ord;
C = col;
end