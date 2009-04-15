function AV = bnActivity(F)

% AV = bnActivity(F) - activity vectors of a Boolean network
% This function uses the truth tables stored as columns in the matrix F of a Boolean
% network and computes the vector of activities for each node (gene).
% The structure of the matrix F is the same as in pbnRnd. The
% activity vectors are stored as rows in the output matrix AV.
% This function uses a fast spectral algorithm and sparse matrix
% structures.
%
% Ilya Shmulevich March 06, 2001 (original version)
% 08/25/05 modified by IS

n = log2(size(F,1));						% maximum number of variables
ng = size(F,2);                             % number of genes
nv = zeros(1,ng);                           % number of inputs for each gene

for i = 1:ng                                % extract number of inputs
    minfind = min(find(F(:,i) == -1));
    if ~isempty(minfind)
        nv(i) = log2(minfind - 1);
    else
        nv(i) = n;
    end
end

AV = zeros(ng,n);
twon = 2.^[0:n];							% we will need powers of 2

for j = 1:ng                                % for each gene
    for i = 1:nv(j),						% this gets a bit ugly
        AV(j,i) = sum(mod(kron(kron(speye(twon(i)),[1 1]),...
            speye(twon(nv(j)-i+1)))*F(1:2^nv(j),j),2))/twon(nv(j))';
    end
end