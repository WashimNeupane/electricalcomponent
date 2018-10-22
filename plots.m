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
spy(RCM);
title("A RCM");

subplot(mm,nn, 4);
spy(banded);
title("Banded A");

subplot(mm,nn, 6);
spy(sprse);
title("Sparse A");

subplot(mm,nn, 5);
spy(AMD);
title("AMD A");

%%plotting temperature gradient ---------------
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
spy(CSR_chol);
title("CSR cholesky");
