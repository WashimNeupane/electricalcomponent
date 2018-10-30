function X = substitution(type,A,B)
if(type == "forward")
X = forward_substitution(A,B);
else
X = backward_substitution(A,B);
end
end

function X = forward_substitution(A, B, divide)
% BACKWARD_SUBSTITUTION Backward substitution
% X = backward_substitution(A, B) solves A*X = B for upper triangular A.
% X = backward_substitution(A, B, false) solves A*X = B for upper triangular A
% without diagonal division.
% If only 2 input arguments, default value for divide is true
if nargin == 2
divide = true;
end
n = size(A,1);
X = B;
for i = 1:n
for j = 1:i-1
X(i) = X(i) - A(i,j) * X(j);
end
if divide
X(i) = X(i) / A(i,i);
end
end
end


function X = backward_substitution(A, B, divide)
% BACKWARD_SUBSTITUTION Backward substitution
% X = backward_substitution(A, B) solves A*X = B for upper triangular A.
% X = backward_substitution(A, B, false) solves A*X = B for upper triangular A
% without diagonal division.
% If only 2 input arguments, default value for divide is true
if nargin == 2
divide = true;
end
n = size(A,1);
X = B;
for i = n:-1:1
for j = i+1:n
X(i,:) = X(i,:) - A(i,j) * X(j,:);
end
if divide
X(i,:) = X(i,:) / A(i,i);
end
end
end
