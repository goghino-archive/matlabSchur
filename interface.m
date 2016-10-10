% Optimal control system
A = '/home/kardos/block_solve/matrices/PardisoMat_10_2.csr';
% Grid size
N = 10;
% Number of scenarios
NS = 2;
%number of worker processes
np = 1;

X = mexsolve(A, N, NS, np);
