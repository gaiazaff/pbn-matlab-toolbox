function [F,varF,nv] = bnThresh2F(M)

% This function constructs a Boolean network in the form of [F,varF,nv]
% (see bnA) from an adjacency matrix M that represents a Boolean threshold
% network. For instance, see:
% PNAS 2004 Apr 6;101(14):4781-6.
%
% The method also includes self-degradation loops so that proteins are
% degraded at the next time point if they have no negative
% (inhibitory) inputs.
% ilya - 07/31/05
% this is a test

n = size(M,1);
if n ~= size(M,2)
    error('M must be square.')
end

% initialization
nv = zeros(1,n);                        % connectivity
for i = 1:n,                            % for each gene i
    nv(i) = nnz(M(:,i));                % number of inputs of gene i
end
nv = nv + 1;                            % add self-input to all genes

F = -1*ones(2^max(nv),n);
varF = -1*ones(max(nv),n);

for i = 1:n,                                            % for each gene i
    varF(1:nv(i),i) = vertcat(i,find(M(:,i)));          % input variables
    if nv(i)>1                                          % if there is at least one non-self input
        S = ff2n(nv(i)-1);                              % all possible inputs (except self)
        s = S*nonzeros(M(:,i));                         % compute weighted sum
        f0 = zeros(size(s));                            % truth table
        f0(s>0) = 1;
        f0(s<0) = 0;
        f0(s==0) = 0;
        f1 = f0;
        f1(s==0) = 1;
        if ~any(M(:,i)<0)                               % if there are no negative inputs
            f1(1) = 0;                                  % self-degradation
        end
        F(1:2^nv(i),i) = vertcat(f0,f1);
        else                                            % if there are no inputs
            F(1:2,i) = [0;0];                           % all-zero truth table
    end
end
