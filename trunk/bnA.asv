function [A,v] = bnA(F,varF,nv,p,method)

% [A,v] = bnA(F,varF,nv,p,method) - state transition matrix of a Boolean network
% This function creates the state transition matrix A of a Boolean network.
% The inputs F, varF and nv are the same as for PBNs (see pbnA) Parameter p
% is a perturbation parameter. It can be used to treat the Boolean network
% as a Markov chain. Another output, v, is the steady-state distribution
% vector of the corresponding Markov chain. The last parameter (optional)
% specifies the method of perturbation; if method == 'cube', then the
% n-cube method is used if method == 'full', then a perturbation p is added
% to all entries of the transition matrix; default is 'full'

% Functions used: bnNextState.m

% Ilya Shmulevich; 04/17/02; 05/09/02
% Modified May 13, 2003 by HL.

n = length(nv); % number of genes
A = zeros(2^n,2^n); % Initialize the output.
b = 2.^[n-1:-1:0]'; % Used to convert binary number to decimal numbers.
bits = [n:-1:1];


if nargin < 5 % Set the default value for the parameter 'method'.
    method = 'full';
end

if lower(method) == 'full'  % depending on what method we use, it will tell us how to scale the matrix A
    z = 2^n;                % if it is 'full', then we add p to all elements and thus, we should subtract 2^n p's from 1.
else
    z = n;                  % if is is the n-cube method, then we only add n p's
end

for i=1:(2^n)
    x = bitget(i-1,bits); % The current row (minus one) as binary number.
    y = bnNextState(x,F,varF,nv);
    j = y*b + 1; % Column index.
    A(i,j) = 1-z*p;
end % for i=1:(2^n)


if lower(method) == 'cube'
    P = [0 1;1 0];      % this is the n-cube method
    for i = 1:n-1,
        P = [P eye(2^i);eye(2^i) P];
    end
else
    P = ones(size(A));  % this is the full method
end

P = P*p;

A = (A + P);

if nargout == 2
   [v,d] = eig(A');
   [maxeig,idxeig] = max(diag(d));		% largest eigenvalue
   v = v(:,idxeig);
   v = v/sum(v);
end
