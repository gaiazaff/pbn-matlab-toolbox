function [bestgene,m] = pbnStateVisit(x,y,A,k0)

% [bestgene,m] = pbnStateVisit(x,y,A,k0) - best gene for intervention
% This function can be used to decide which gene is the best possible
% target for intervention if we want to go from state x to state y. Vectors
% x and y are binary vectors (1 by n) representing starting and destination
% states. Matrix A is the state transition matrix as produced by pbnA (see
% description) k0 is a user-selectable parameter; we are interested in the
% probability that the network, starting in x will visit state y before
% time k0. The output, bestgene, is the index of the best candidate gene
% with which to intervene. In other words, if we flip that gene, we are
% most likely to get to our destination y before time k0. Another optional
% output is m, which is the mean first passage times. The size of m is 1xn
% and each element corresponds to the mean first passage time corresponding
% to the perturbation of each gene. The function also plots the
% probabilities for every k = 1...k0, for every possible gene (1,...,n).

% Ilya Shmulevich; Aug 27, 2001; Nov 01, 2001

n = length(x);											% number of genes
f = zeros(2^n,k0);									% this will store F_k(x,y)

xdec = b2d(x) + 1;									% convert x and y to decimal
ydec = b2d(y) + 1;									% and add 1 for indexing

f(:,1) = A(:,ydec);									% for k = 1, it's just the column of A corresp. to y
Q = A;													% Q is A with that column
Q(:,ydec) = zeros(2^n,1);							% replaced by zeros

for k = 2:k0,											% now use recursive formula
   f(:,k) = Q*f(:,k-1);
end

pertgene = zeros(1,n);								% this will contain the n possible starting states

for i = 1:n
   tempstate = x;
   tempstate(i) = ~tempstate(i);					% flip the i-th bit
   pertgene(i) = b2d(tempstate)+1;				% decimal version of that state
end

H = cumsum(f(pertgene,:)');						% sum F_k(.,y) from 1..k0
plot(H)													% and plot the result
legend(num2str([1:n]'),0)

[tmp,bestgene] = max(H(k0,:));					% output best gene candidate for intervention

m = sum((repmat([1:k0],n,1) .* f(pertgene,:))');			% mean first passage times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = b2d(x)
% Y = BIN2DEC (X)
% This function can be used to convert a binary vector to a decimal integer.

p = 2.^[length(x)-1:-1:0];
y = x * p';