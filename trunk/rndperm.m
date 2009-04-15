function y = rndperm(n,low,high)

% y = rndperm(n,low,high) - random permutation of a set of integers
% rndperm(n,low,high) returns a random permutation of integers between low
% and high (inclusing limits). Restrictions: low<=high and 0<=n<=high-low+1.

% Functions used:

% 14.04.2003 by HL.

% Check the restrictions.
if low>high
   error('low must satisfy low<=high.');
end % if low>high
if n>(high-low+1) | n<0
   error('n must satisfy 0<=n<=high-low+1.');
end % if n>(high-low+1) | n<0

% Number of integers between (and including) low and high.
N = high-low+1;

[temp,y] = sort(rand(1,n));
y = y + (low - 1);
