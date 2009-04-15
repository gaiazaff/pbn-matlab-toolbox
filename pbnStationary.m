function pmf = pbnStationary(Y,Yp,p)

% pmf = pbnStationary(Y,Yp,p) - empirical stationary distribution
% This function takes the history Y produced either by pbnRun or pbnTemper
% and produces an empirical stationary distribution, by counting occurences
% of each state There are two more optional input variables: Yp and p Yp is
% produced along with Y by pbnTemper() and contains the history of the
% temperatures p is the desired temperature to extract from Y. Example:
% pbnStationary(Y,Yp,0.2) will extract only those states from Y which
% correspond to temperature 0.2

% Ilya Shmulevich; Aug. 22, 2001

n = size(Y,2);										% number of genes

if nargin == 3
   Y = Y(find(Yp == p),:);						% extract those states corresponding to p
end

p = 2.^[n-1:-1:0];								% convert to decimal
pmf = hist(Y * p',2^n);							% and create histogram
pmf = pmf / sum(pmf);							% then normalize