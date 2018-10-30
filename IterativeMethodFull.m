function x1 = IterativeMethod(IterationType,a,b,omega)
n = length(b);
x_ny = zeros(1,n);
x0 = zeros(n,1);
tol=1;
max = 100;

if nargin<4
    omega = 1;
end

%%Jacobi Iteration
if(IterationType == "Jacobi")
    x = zeros(n,1);
for i = 1:n 
     x(i) = ((b(i) - a(i,[1:i-1,i+1:n]) * x0([1:i-1,i+1:n])) / a(i,i));
end
x1 = x;
k = 1;

while norm(x1-x0,1) > tol
    for j = 1 : n
     x_ny(j) = ((b(j) - a(j,[1:j-1,j+1:n]) * x1([1:j-1,j+1:n])) / a(j,j));
    end
    x0 = x1;
    x1 = x_ny';
    k = k + 1;
end

end

%Gauss_seidel iteration

if(IterationType == "GaussSeidel")
x1 = SOR(1);
end

%SOR iteration
if(IterationType == "SOR")
x1 = SOR(omega);
end


%Conjugate gradient method
if(IterationType == "ConjugateGradient")

%for first iteration
x = x0;
r = b - (a* x );  
d = -r;
Ad = a*d;
dtAd = d'*Ad;
mag = (r'*d)/dtAd;
x = x + mag*d;

while(norm(r)>tol)  
    for k = 1:n
       r = r - mag*Ad;       
       B = (r'*Ad)/dtAd;
       d = -r + B*d;
       Ad = a*d;
       dtAd = d'*Ad;
       mag = (r'*d)/dtAd;
       x = x + mag*d;
    end
end
 x1 = x;
end
end

%function SOR
function x1 = SOR(omega)
n = length(b);
x_ny = zeros(1,n);
x0 = zeros(n,1);
tol=1;
max = 100;

    x = zeros(n,1);
        for j = 1:n 
             if j==1
                 x(j) = ((b(j) - a(j,[1:j-1,j+1:n]) * x0([1:j-1,j+1:n])) / a(j,j));
             else    
                  x(j) = ((1-omega)*x(j-1))+ ((b(j) - a(j,[1:j-1,j+1:n]) * x0([1:j-1,j+1:n]))* (omega / a(j,j)));     
             end
                  x_pre(j) = x0(j);
                  x0(j) = x(j);
        end
            x1 = x;
            k = 1;
            
        while k<max || norm(x1-x_pre,1) > tol
             for j = 1 : n
                if j == 1
                     x_ny(j) = ((b(j) - a(j,[1:j-1,j+1:n]) * x0([1:j-1,j+1:n])) / a(j,j));
                else
                     x_ny(j) = ((1-omega)*x_ny(j-1))+((b(j) - a(j,[1:j-1,j+1:n]) * x1([1:j-1,j+1:n]))* (omega / a(j,j)));             
                end
                x_pre(j) = x1(j);
                 x1(j) = x_ny(j);
             end
             x0 = x1;
             x1 = x_ny';
             k = k + 1;
        end
end



