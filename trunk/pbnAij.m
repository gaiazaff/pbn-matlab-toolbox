function pij = pbnAij(F,varF,nf,nv,cij,i,j,p)

% pij = pbnAij(F,varF,nf,nv,cij,i,j,p) - transition prob. between two states
% This function returns the ij-th entry (probability) in the state
% transition matrix A corresponding to a Probabilistic Boolean Network
% (PBN). The output pij is a probability.
% INPUT:
% i,j   - The index of interest into the A matrix (starting from 1)
% p     - The probability of a gene being randomly perturbed
% Other inputs are defined e.g. in pbnRnd.m

% Ilya Shmulevich; Aug. 3, 2002
% Modified May 13, 2003 by HL.

n = length(nf); % number of genes
cnf = [0 cumsum(nf)];

j = j - 1; % indexing starts from 1
i = i - 1;
ibin = bitget(i,n:-1:1); % i and j as binary numbers.
jbin = bitget(j,n:-1:1);
pij = 1; % initialize pij
b = 2.^[size(varF,1)-1:-1:0]';

for g=1:n % for each target gene
    cijtotal = 0;
    for m=1:nf(g)
        c = cnf(g)+m; % Column index.
        % The probability that the current function will output the desired
        % output value for the current node.
        y = F(ibin(varF(1:nv(c),c))*b(end-nv(c)+1:end)+1,c);
        cijtotal = cijtotal + cij(m,g)*(y == jbin(g));
    end
    pij = pij * cijtotal;
end

if p > 0 % if some genes can be perturbed
    hamdist = sum(bitxor(ibin,jbin));
    pij = (pij * (1-p)^n) + p^hamdist * (1-p)^(n-hamdist) * (i~=j);
end
