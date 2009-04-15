function [P,FA] = bnAttractorProximity(ab)
% SYNTAX: P = bnAttractorProximity(ab)
% This function creates a proximity matrix P between the attractors. The
% only input is the vector ab, produced by bnAttractor. The distance
% between attractor i and j is stored in P(i,j) and is computed as follows:
% For every attractor, we determine a vector of length equal to the number
% of genes which contains the fraction of time that each gene is ON in that
% attractor. Then, the distance between two attractors is the Euclidean
% distance between the corresponding vectors.
% An optional output parameter, FA, will output the real-valued vectors (as
% rows) corresponding to each attractor.
% Ilya Shmulevich; 12/13/02

numattractors = length(unique(ab(ab<0)));               % number of attractors
n = log2(length(ab));                                   % number of genes
FA = zeros(numattractors,n);                            % each row will correspond to an attractor and contain the fraction of time that each gene is 1
P = zeros(numattractors,numattractors);                 % matrix of proximities


for attr = 1:numattractors,                             % for each attractor
    FA(attr,:) = mean(dec2binarr(find(ab == -attr) - 1,n),1);   % compute the average fractional "gene activity" vector
end

for i = 1:numattractors,
    for j = i:numattractors,
        P(i,j) = norm(FA(i,:)-FA(j,:));                 % compute distances between rows of FA
    end
end
P = P + tril(P',-1);                                    % since we only compute half of the matrix, make it symmetric

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function will convert a decimal to a binary array
% e.g. dec2binarr([0:7],3) produces all binary combinations for 3 variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binarr=dec2binarr(dec,len)

binarr=zeros(length(dec),len);

for i=1:length(dec)
    binarr(i,:) = bitget(dec(i),len:-1:1);
end  