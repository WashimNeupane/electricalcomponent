function [x1,numberOfIterations,runtime] = IterativeMethodCSR(type, a,row,col,b,omega)
tic
if nargin < 6
    omega =1;
end

nb = length(b);
n = length(a);
nr = length(row);
x0 = zeros(nr-1,1);
numberOfIterations =1;
x = zeros(nb,1);
tol =1*10;

%Jacobi
if(type == "Jacobi")
 for i =1:nr-1         
        bb = b(i);
        for j = row(i):row(i+1)-1
        if(col(j) == i)
            ii = j;
        end
        end
        aa = a(1,[row(i):ii-1, ii+1:row(i+1)-1]);        
        xx =  col([row(i):ii-1, ii+1:row(i+1)-1])';
        xx = x0(xx,1);
                
        x(i) = (bb - aa*xx)/a(ii);        
    end
    
x1 = x;
while norm(x1-x0,1)>tol
    for i = 1 : nr-1        
        bb = b(i);        
        for j = row(i):row(i+1)-1
        if(col(j) == i)
            ii = j;
        end
        end
               
        aa = a(1,[row(i):ii-1, ii+1:row(i+1)-1]);
        xx =  col([row(i):ii-1, ii+1:row(i+1)-1])';
        xx = x1(xx,1);                
        x_ny(i) = (bb - aa*xx)/a(ii); 
    end
    x0 = x1;
    x1 = x_ny';
    numberOfIterations = numberOfIterations + 1;
end

elseif(type == "GaussSiedel")
[x1,numberOfIterations] = SOR(a,row,col,b,1,tol);

elseif(type == "SOR")
[x1,numberOfIterations] = SOR(a,row,col,b,omega,tol);

elseif(type == "ConjugateGradient")
n = length(a);
nr = length(row);
x0 = zeros(nr-1,1);
tol =1*10^-10;

%for first iteration
x = x0;
for i=1:nr-1 
   cc = (a(1,row(i):row(i+1)-1));
   cc2 = (col(1,row(i):row(i+1)-1))';
   xx = x(cc2,1);
   r = b - cc * xx; 
end
 
d = -r;
Ad = zeros(3,1);

for i=1:nr-1
   cc = a(1,row(i):row(i+1)-1);
   cc2 = (col(1,row(i):row(i+1)-1))';
   xx = d(cc2,1);
   Ad(i) = cc*xx; 
end

dtAd = d'*Ad;
mag = (r'*d)/dtAd;
x = x + mag*d;

while(norm(r)>tol)  
    for k = 1:n
       r = r - mag*Ad;       
       B = (r'*Ad)/dtAd;
       d = -r + B*d;
       
       for i=1:nr-1
        cc = a(1,row(i):row(i+1)-1);
        cc2 = (col(1,row(i):row(i+1)-1))';
        xx = d(cc2,1);
        Ad(i) = cc*xx; 
       end
       
       dtAd = d'*Ad;
       mag = (r'*d)/dtAd;
       x = x + mag*d;
    end
    numberOfIterations = numberOfIterations + 1;
end
 x1 = x;
end
runtime = toc;
end


function [x1,k] = SOR(a,row,col,b,omega, tol)
nb = length(b);
n = length(a);
nr = length(row);
x0 = zeros(nr-1,1);
k =1;
x = zeros(nb,1);

 for i =1:nr-1
        bb = b(i);
        for j = row(i):row(i+1)-1
        if(col(j) == i)
            ii = j;
        end
        end
        aa = a(1,[row(i):ii-1, ii+1:row(i+1)-1]);        
        xx =  col([row(i):ii-1, ii+1:row(i+1)-1])';
        xx = x0(xx,1); 
        
       if (i==1)
            x(i) = (bb - aa*xx)/a(ii);
       else
           x(i)= ((1-omega)*x(i-1)) +  (bb - aa*xx)* (omega/a(ii)); 
       end
       x_pre(i) = x0(i);
       x0(i) = x(i);
      end
x1 = x;

while norm(x1-x_pre',1)>tol
     for i =1:nr-1
        bb = b(i);
        for j = row(i):row(i+1)-1
        if(col(j) == i)
            ii = j;
        end
        end
        aa = a(1,[row(i):ii-1, ii+1:row(i+1)-1]);        
        xx =  col([row(i):ii-1, ii+1:row(i+1)-1])';
        xx = x1(xx,1); 
        
       if (i==1)
            x_ny(i) = (bb - aa*xx)/a(ii);
       else
           x_ny(i)= ((1-omega)*x_ny(i-1)) +  (bb - aa*xx)* (omega/a(ii)); 
       end
       x_pre(i) = x0(i);
       x1(i) = x_ny(i);
      end
    x0 = x1;
    x1 = x_ny';
    k = k+1;
end
end
