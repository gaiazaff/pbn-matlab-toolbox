function P = bnAttractorMarkov(ab,F,varF,nv,sm)

% P = bnAttractorMarkov(ab,F,varF,nv,sm) - attractor proximity matrix
% This function computes an attractor proximity matrix, P, for a Boolean
% network. For every state in an attractor, we transiently perturb every
% gene and determine to which attractor (including possibly itself) the 
% network will flow. We tabulate the number of times this will occur and
% P(i,j)  contains the total number of times the network will transition, 
% when perturbed, from attractor i to attractor j. Thus, P is a square 
% matrix with the number of rows/columns equal to the number of attractors.
% The inputs to the functions are ab (produced by bnAttractor) and F, varF
% and nv (which define the Boolean network). Note that in order to convert
% P to a real stochastic matrix, we can pass a fifth, optional argument
% that can be equal to anything.

% Functions used: bnNextState.m

% Ilya Shmulevich; 12/11/02
% Modified May 14, 2003 by HL.

numattractors = length(unique(ab(ab<0)));               % number of attractors
n = log2(length(ab));                                   % number of genes
twon = 2.^[n-1:-1:0]';
P = zeros(numattractors,numattractors);                 % "transition" matrix - see Fig 12.12 (p. 490) in Origins of Order

for attr = 1:numattractors,                             % for each attractor
    attractor_states = find(ab == -attr) - 1;           % states of the current attractor
    for i = 1:length(attractor_states),                 % for each state in that attractor
        for h = 1:n,                                    % perturb each gene
            jumpto = ab(bnNextState(bitget(bitxor(attractor_states(i),twon(h)),n:-1:1),F,varF,nv)*twon+1);    % attractor/basin to which we jump
            P(attr,abs(jumpto)) = P(attr,abs(jumpto)) + 1;  % abs() is because attractors are negative and basins are positive
        end
    end
end

if nargin == 5                                          % if we want a stochastic matrix
    P = P./repmat(sum(P')',1,numattractors);
end