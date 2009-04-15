function D = DistortionVectors(d,e,i,D,c,emax)

% D = DistortionVectors(d,e,i,D,c,emax) - Returns proper distortion vectors
%
% DistortionVectors -function returns proper distortion vectors in D (or their
% permutations) such that together with the optimal function, or its truth
% table, they define the set of functions {f : e(f)<=emax} that have error
% size smaller than equal to (ref. Best-Fit Extension Problem) emax. See a
% paper for details.
%
% Input variables:
% d     - a (current permuted) distortion vector
% e     - a (current) error size
% i     - index of the right most bit that is set in the previous (upper) level
%         the recursive algorithm
% D     - set of proper (permuted) distortion vector. Each row of D corresponds
%         to a single distortion vector.
% c     - (permuted) "weight" vector
% emax  - maximum error size
%
% Output variables:
% D     - set of proper (permuted) distortion vector. See above for details.
%
% Functions used:

% 08.04.2002 by Harri Lähdesmäki
% Modified:

j = i + 1;
stop = length(d); % = 2^k
while j<=stop & (e+c(j))<=emax
    d(j) = 1;
    D = [D;d];
    D = DistortionVectors(d,e+c(j),j,D,c,emax);
    d(j) = 0;
    j = j+1;
end % while j<=stop & (e+c(j))<=emax
