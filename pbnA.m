function [A,v] = pbnA(F,varF,nf,nv,cij,p)

% [A,v] = pbnA(F,varF,nf,nv,cij,p) - state transition matrix of a PBN
% This function creates the 2^n x 2^n transition matrix A corresponding to
% a PBN. Another optional output argument of the function is the stationary
% distribution v.
% INPUT:
% p     - probability of a gene being randomly perturbed
% Other inputs are defined e.g. in pbnRnd.m

% Functions used: pbnAij.m

% Ilya Shmulevich; Aug. 14, 2001
% Modified May 13, 2003 by HL.

n = length(nf); % number of genes
A = zeros(2^n,2^n);

% Set each element of A.
for i = 1:2^n,
    for j = 1:2^n,
        A(i,j) = pbnAij(F,varF,nf,nv,cij,i,j,p);
    end
end

if nargout == 2
    [v,d] = eig(A');
    [maxeig,idxeig] = max(diag(d));		% largest eigenvalue
    v = v(:,idxeig);
    v = v/sum(v);
end
