function y = pmfrnd(pmf)

% y = pmfrnd(pmf) - sample a discrete distribution
% This function takes as input a vector pmf, which is a probability mass
% function (table) e.g. if pmf = [.3 .2 .1 .4] it outputs an index (e.g. 1,2,3
% or 4) with the corresponding probabilities.
% It is *assumed* that sum(pmf)==1

% Ilya Shmulevich; Aug. 14, 2001

cdf = cumsum(pmf);			% make cumulative distr. function
u = rand;						% pick a random number in [0,1]
y = sum(u > cdf) + 1;		% find where it falls