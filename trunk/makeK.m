function K = makeK(Li)
% K = makeK(Li) - matrix K
% This function creates the matrix K (see PBN paper) It is of size N x n Li
% is a vector [l(1) l(2) .... l(n)], where l(i) is the number of possible
% functions for gene i

% Ilya Shmulevich Aug. 16, 2001
% Modified May 14, 2003 by HL.

N = prod(Li);							% number of network realizations
n = length(Li);						% number of genes
K = zeros(N,n);						% initialize K
idx = 1;									% index into K
[y,K,idx] = makeK3(ones(1,n),Li,1,K,idx);				% call a recursive function

%---------------------------
function [y,K,idx] = makeK3(x,Li,ptr,K,idx)

n = length(Li);						% number of genes

for i = 1:Li(ptr),
   x(ptr) = i;
   if (ptr == n)
      K(idx,:) = x;
      idx = idx + 1;
   end
   if ptr < n
      [x,K,idx] = makeK3(x,Li,ptr+1,K,idx);
   end
end

y = x;