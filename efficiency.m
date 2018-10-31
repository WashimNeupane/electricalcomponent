%%test matrices
tolerance = [1, 0.1, 0.001, 0.00001, 1*10^-10, 1* 10^-15];

%%workspaceEfficiency.mat contains a matrix of runtimes and iteration
%%each column in runtime represents Jacobi,GaussSeidel,ConjugateGradient,
%%SOR with omega 0.1,0.5,0.94,0.77 respectively.
load('workspaceEfficiency.mat')

%%Gives the memory size associated with each of these matrices.
whos full packed banded sprse CSR CSR_RRow CSR_RCol

%%Gives the memory size of all cholesky factorised matrices
whos full_chol packed_chol_linear banded_chol sprse_chol

%--------------------------------------------------------------------
%%Runtimes times for convergence with changing values of tolerance
figure('rend','painters','pos',[1100 70 700 900]);
subplot(2,1,1)
plot(iteration(1:5,1:7));
title('Number of Iterations for convergence with varying toleance')
xticks([1 2 3 4 5]);
xticklabels({'1*10^-10','1*10^-8','1*10^-4','1*10^-1','10'})
xlabel('tolerance'); ylabel('number of iterations')
legend('Jacobi','Gauss Siedel','Conjugate Gradient','SOR-Omega = 0.1','SOR-Omega = 0.5','SOR-Omega = 0.94','SOR-Omega = 0.77');

subplot(2,1,2)
plot(runtime(7:11,1:7));
title('Runtimes for convergence with varying toleance')
xticks([1 2 3 4 5]);
xticklabels({'1*10^(-10)','1*10^(-8)','1*10^(-4)','1*10^(-1)','10'})
xlabel('tolerance'); ylabel('time in sec')
legend('Jacobi','Gauss Siedel','Conjugate Gradient','SOR-Omega = 0.1','SOR-Omega = 0.5','SOR-Omega = 0.94','SOR-Omega = 0.77');
%--------------------------------------------------------------------

%%Plot runtime for Jacobi, GaussSeidel and Conjugate Gradient,SOR
%%constant tolerance level
figure('rend','painters','pos',[1100 70 700 900]);
plot(runtime(1:5,1:7));
title('Runtimes for Iterative Methods with constant tolerance')
xticks([1 2 3 4 5 6]);
xticklabels({'1','2','3','4','5','6'})
xlabel('trial number'); ylabel('time in sec')
legend('Jacobi','Gauss Siedel','Conjugate Gradient','SOR-Omega = 0.1','SOR-Omega = 0.5','SOR-Omega = 0.94','SOR-Omega = 0.77');
%----------------------------------------------------------------------

%%Overall times for convergence with changing values of omega for SOR
figure('rend','painters','pos',[1100 70 700 900]);
subplot(2,1,1)
plot(runtime(1:5,4:7));
title('Runtimes for SOR with different omega values with constant tolerance')
xticks([1 2 3 4 5 6]);
xticklabels({'1','2','3','4','5','6'})
xlabel('trial number'); ylabel('time in sec')
legend('SOR-Omega = 0.1','SOR-Omega = 0.5','SOR-Omega = 0.94','SOR-Omega = 0.77');

subplot(2,1,2)
plot(iteration(1:5,4:7));
title('Number of Iterations for SOR with different omega values with varying toleance')
xticks([1 2 3 4 5]);
xticklabels({'1*10^-10','1*10^-8','1*10^-4','1*10^-1','10'})
xlabel('tolerance'); ylabel('number of iterations')
legend('SOR-Omega = 0.1','SOR-Omega = 0.5','SOR-Omega = 0.94','SOR-Omega = 0.77');






