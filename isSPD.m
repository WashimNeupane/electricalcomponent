function bool = isSPD(A)
dim = size(A);
if ~(dim(1)==dim(2))
    error('A must be square')    
end
if (all((all(A)~=all(A'))))
    error('A must be symmetric')
end
if ~(all(eig(A))> 0)
    error('Matrix is at least not positive definite')
end
    bool = true;
end