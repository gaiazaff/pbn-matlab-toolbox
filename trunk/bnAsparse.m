function [A,Avec] = bnAsparse(F,varF,nv)

% [A,Avec] = bnAsparse(F,varF,nv) - sparse state transition matrix of a Boolean network
% This function creates the state transition matrix A of a Boolean network,
% stored as a sparse matrix. An optional output is A in vector form (Avec)
% The inputs F, varF and nv are the same as for PBNs (see pbnA)

% Functions used: bnNextState.m

% Ilya Shmulevich; 9/9/02
% Modified May 13, 2003 by HL.

n = length(nv);         % number of genes

j = zeros(2^n,1);
bits = [n:-1:1];
b = 2.^(bits'-1);          % Used to convert binary number to decimal numbers.

for i = 1:2^n
    x = bitget(i-1,bits); % The state (minus one) as binary number.
    y = bnNextState(x,F,varF,nv);
    d = y*b + 1;        % Output as a decimal number.
    j(i) = d;
end % for i = 1:2^n

i = [1:2^n]';

A = sparse(i,j,1,2^n,2^n);
if nargout == 2
    Avec = j;
end
