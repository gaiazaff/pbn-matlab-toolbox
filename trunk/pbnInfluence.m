function G = pbnInfluence(F,varF,nf,nv,cij,d)

% G = pbnInfluence(F,varF,nf,nv,cij,d) - Influence matrix of a PBN
% This function computes the influences matrix Gamma. G(i,j) contains the
% influence of gene i on gene j.
% INPUT:
% F, varF, nf, nv, cij defines an independent PBN (see e.g. pbnRnd.m)
% d - joint distribution of the variables.
% OUTPUT:
% G - Influence matrix.

% May 14, 2003 by Harri Lähdesmäki

n = length(nf); % The number of nodes.

G = zeros(n,n); % initialize influence matrix

cnf = [0,cumsum(nf)];

for i=1:n % Influencing variable/node.
    
    for j=1:n % Influenced variable/node.
        
        % Extract the parameters for the j:th node.
        Fj = F(:,(cnf(j)+1):(cnf(j)+nf(j))); % Function set (may include -1s).
        varFj = varF(:,(cnf(j)+1):(cnf(j)+nf(j))); % Sets of variables (may include -1s).
        nvj = nv((cnf(j)+1):(cnf(j)+nf(j))); % Number of variables.
        cijj = cij(1:nf(j),j); % Selection probabilities.
        
        % Check all the predictor functions of the j:th node.
        for k=1:nf(j)
            
            vars = varFj(1:nvj(k),k); % Variables of the k:th pred. function of the j:th node.
            f = Fj(1:2^nvj(k),k); % k:th predictor function of the j:th node.
            isthereany = find(vars==i);
            
            % If the i:th node is among the variables of the current function.
            if ~isempty(isthereany)
                
                if length(vars)==1
                    G(i,j) = G(i,j) + (f(1)~=f(2))*cijj(k); % If the only variable is essential.
                else
                    varsminusi = vars;
                    varsminusi(isthereany) = []; % Remove the element having value i, set of variables minus i.
                    %disp([i,varsminusi']);,pause;
                    ff = permutevars(f,vars,[i,varsminusi']); % Change the ordering of the variables, i:th node first.
                    md = margpdf(d,varsminusi); % Joint marginal distribution over the "set of variables minus i."
                    infl = sum((ff(1:2^(nvj(k)-1))~=ff(2^(nvj(k)-1)+1:2^nvj(k)))'.*md); % Influence of the current function.
                    G(i,j) = G(i,j) + infl*cijj(k);
                end % if length(vars)==1
            end % if ~isempty(isthereany)
        end % for k=1:nf(j)
    end % for j=1:n
end % for i=1:n
