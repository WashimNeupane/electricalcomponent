%%plotting all images ---------------
mm = 4; nn=2;
figure('rend','painters','pos',[1100 70 700 900]);
subplot(mm,nn, 1);
spy(A);
title("A");

subplot(mm,nn, 2);
spy(sparse(A));
title("CSR");

subplot(mm,nn, 3);
spy(banded_RCM);
title("A RCM");

subplot(mm,nn, 4);
spy(banded);
title("Banded A");

subplot(mm,nn, 6);
spy(sprse);
title("Sparse A");

subplot(mm,nn, 5);
spy(sparse_AMD);
title("AMD A");

%PLOT TEMPERATURE MAP--------------------------------------------
[tri,temperatures,x,y] = makegrid(T);
figure('rend','painters','pos',[1100 70 700 900]);
trisurf(tri,x,y,temperatures)
title("temperature plot");
view(2)
shading interp
colormap hot
colorbar
axis image
grid off


%%plotting cholesky ---------------
mm = 2; nn=2;
figure('rend','painters','pos',[1100 70 700 900]);
subplot(mm,nn, 1);
spy(packed_chol);
title("packed cholesky");

subplot(mm,nn, 2);
spy(banded_chol);
title("banded cholesky");

subplot(mm,nn, 3);
spy(sprse_chol);
title("sparse cholesky");

subplot(mm,nn, 4);
% spy(CSR_chol);
title("CSR cholesky");


%%Functions
function [tri,temperatures,x,y]= makegrid(T)
[X,Y]= meshgrid(0:.01:0.06,0.06:-.01:0);
temperatures = T(~isnan(T));
x = X(~isnan(T));
y = Y(~isnan(T));
tri = delaunay(x, y);

for ii = 1:size(tri,1)
indices = tri(ii,:);
vertices = [x(indices) y(indices)];
centroid = mean(vertices);
text(centroid(1),centroid(2), num2str(ii));
end
tri([29,44,36],:) = [];
end
