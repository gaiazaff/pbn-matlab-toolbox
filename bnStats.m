function s = bnStats(ab,d)

% s = bnStats(ab,d) - statistics of a Boolean Network
% This function takes as input vectors ab and d, as produced by
% bnAttractor, and outputs a structure s.
% The structure contains:
% s.numattractors:      the number of attractors in the network
% s.attractor_sizes:    the sizes of each attractor
% s.basin_sizes:        the sizes of the corresponding basins
% s.delta_k:            the maximum distance in each basin (to the
%                       attractor)
% s.size_Dj:            cell array: the jth element in the kth cell
%                       contains the number of states in the kth basin
%                       whose distance is <= j

% Ilya Shmulevich; Harri Lähdesmäki; 09/12/02
% Modified by Ilya (04/29/04): fixed bug related to empty basins
% Modified by Harri (05/05/06): fixed bug related to empty basins

s.numattractors = length(unique(ab(ab<0)));                 % number of attractors
[n,x] = hist(ab(find(ab<0)),s.numattractors);
s.attractor_sizes = fliplr(n);
if length(find(ab>0))>0
%  [n,x] = hist(ab(find(ab>0)),s.numattractors);
  [n,x] = histc(ab(find(ab>0)),1:s.numattractors);
else
  n = zeros(size(s.attractor_sizes));
end
s.basin_sizes = n;
s.delta_k = zeros(1,s.numattractors);
s.size_Dj = cell(1,s.numattractors);

for k = 1:s.numattractors,                                  % this creates the delta_k vector
    if ~isempty(find(ab == k))                              % if there is a basin
        sdeltak = max(d(find(ab == k)));
    else
        sdeltak = 0;                                        % in case there is no basin
    end
    s.delta_k(k) = sdeltak;
    Dkj = length(find(ab == -k))*ones(s.delta_k(k)+1,1);
    for j = 1:s.delta_k(k),
        Dkj(j+1) = Dkj(j+1) + length(find(d(find(ab == k))<=j));
    end
    s.size_Dj{k} = Dkj;
end
