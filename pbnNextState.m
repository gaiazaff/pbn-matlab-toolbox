function y = pbnNextState(x,F,varF,nf,nv,cij,p)

% y = pbnNextState(x,F,varF,nf,nv,cij,p) - one step of a PBN
% This function performs one step of the PBN. The current state is x
% (binary vector of 1 x n, where n is # of genes) For other parameters, see
% e.g. pbnrnd.m

% Ilya Shmulevich; Aug. 14, 2001
% Modified May 12, 2003 by HL.

n = length(x);          % number of genes
y = zeros(size(x));     % initialize next state
gam = rand(1,n) < p;    % generate a random perturbation vector

if any(gam) % if gam is not all zeros
    y = bitxor(x,gam); % perturb some genes
else
    cnf = [0,cumsum(nf)];
    b = 2.^[size(varF,1)-1:-1:0]';
    
    % Update each node separately.
    for i=1:n
        pmf = cij(1:nf(i),i); % extract its cij probabilities
        j = pmfrnd(pmf); % pick a predictor at random (note: 1 <= j <= nf(i))
        k = cnf(i)+j; % Index of the random selected predictor.
        y(i) = F(x(varF(1:nv(k),k))*b(end-nv(k)+1:end)+1,k); % Output value for the selected function.
    end % for i=1:n
end % if any(gam)
