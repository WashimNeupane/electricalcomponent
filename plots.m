%%plotting all images ---------------
figure('rend','painters','pos',[1100 70 700 900]);
subplot(2,2, [1 2]);
spy(A);
title("A");
 
subplot(2,2, 3);
spy(banded_RCM);
title("Banded A");

subplot(2,2, 4);
spy(sprse);
title("Sparse A");


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
figure('rend','painters','pos',[1100 70 700 900]);
subplot(2,2, 1);
spy(full_chol);
title("full cholesky");

subplot(2,2, 2);
spy(packed_chol);
title("packed cholesky");

subplot(2,2, 3);
spy(banded_chol_unpacked);
title("banded cholesky");

subplot(2,2, 4);
spy(sprse_chol);
title("sparse cholesky");

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
