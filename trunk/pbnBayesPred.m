function [Fhat,varF,errF] = pbnBayesPred(v,k,F)

% [Fhat,varF,errF] = pbnBayesPred(v,k,F) - Bayes predictors for a PBN
%
% Function finds the best (Bayes) steady-state predictors for each node in
% the network (the node itself is not allowed to be used as predictor node
% for itself). In case of tie, a random selected function with minimum
% error-size is returned. 
% INPUT:
% v - The stationary distribution of a PBN, i.e., f(1) = Pr{X = 00...00},
%     f(2) = Pr{X = 00...01}, ..., f(2^n) = Pr{X = 11...11}.
% k - The number of variables used in each predictor function.
% F - The set of Boolean predictors to be used in the inference. F is a
%     (2^k)-by-nf binary matrix, where k is the number of variables in each
%     function and nf is the number of functions. If F is an empty matrix,
%     then the function class is considered to contain all k-variable
%     Boolean functions. Let f = F(:,j) be the j:th column of F. Then, f(0)
%     defines the output value for input vector 00...00, f(1) for 00...01,
%     f(2) for 00...10, ..., and f(2^k) for 11...11. Input vectors are
%     interpreted such that the k:th (left most) bit defines the value for
%     the first input node/variable, the k-1:th bit defines the value for
%     the second input node, ..., and the first (right most) bit defines
%     the value for the last input node
% OUTPUT:
% Fhat  - A matrix of the Bayes predictors for the all nodes. Fhat has size
%         (2^k)-by-n, where n is the number of nodes in the PBN. Fhat(:,i)
%         defines the best functions for the i:th node. Each column in Fhat
%         is interpreted as the columns in F (see above).
% varF  - A k-by-n matrix which defines the predictor variables of the best
%         predictors. varF(:,i) defines the predictors for the i:th
%         function in F.
% errF  - The true error of the Bayes predictors for all the target nodes.
%         errF is a 1-by-n vector.

% Functions used: margpdf.m randintex.m

% 06.05.2003 by Harri Lähdesmäki


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define and initialize some variables.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

n = round(log2(length(v))); % The number of nodes.
kk = 2^k; % Needed often.

combnum = nchoosek(n-1,k); % The number of different variable combinations.

% All variable combinations will be generated in advance. Check that he
% number of combinations is "samll" enough. This will work only for small
% enough PBNs.
if combnum>20000 % Limit the number of possible combinations.
    error('Too many variable combinations...')
end % if combnum>20000

starti = 1;
stopi = combnum;

% Initialize the output matrices.
Fhat = zeros(kk,n);
varF = zeros(k,n);
errF = zeros(1,n);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The main loop separately for unconstrained and constrained case.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% If F is an empty matrix, then the function class is considered to contain
% all k-variables Boolean functions (unconstrained case).
if isempty(F)
    
    % Run through all the nodes.
    for i=1:n
        
        minerr = realmax; % Keep track of the minimum error.
        
        % Generate all variable combinations in advance. Remove the current
        % node (target node) from the set of predictors.
        IAll = nchoosek([1:i-1,i+1:n],k);
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Run through all variable combinations.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for j=starti:stopi
            
            % The current variable combinations (in lexicographical ordering).
            I = IAll(j,:);
            
            % Compute the corresponding marginal distribution (including
            % the target node itself).
            vv = margpdf(v,[i,I]);
            
            % Find the optimal (Bayes) k-variable predictor for the current
            % input variable combination.
            [OptErr,OptF] = min([vv(1:kk);vv(kk+1:2*kk)]);
            OptF = 2 - OptF;
            % All output bits having tie are set uniformly randomly.
            Ties = (vv(:,1:kk)==vv(:,kk+1:2*kk));
            OptF(Ties) = (rand(1,sum(Ties(:))))>0.5;
            % The corresponding probability of error.
            OptErr = sum(OptErr);
            
            % Keep track of the minimum error.
            if OptErr<minerr
                Fhat(:,i) = OptF';
                varF(:,i) = I';
                errF(i) = OptErr;
                minerr = OptErr;
            end % if OptErr<minerr
            
        end % for i=starti:stopi
    end % for i=1:n
    
else
    
    % If F is non-empty matrix, then the function class is defined by F.
    
    nf = size(F,2); % The number of functions in F.
    Oneskknf = ones(kk,nf);
    Zeroskknf = zeros(kk,nf);
    Ones1nf = ones(1,nf);
    
    % Run through all the nodes.
    for i=1:n
        
        minerr = realmax; % Keep track of the minimum error.
        
        % Generate all variable combinations in advance. Remove the current
        % node from the set of predictors.
        IAll = nchoosek([1:i-1,i+1:n],k);
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Run through all variable combinations.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for j=starti:stopi
            
            % The current variable combinations (in lexicographical ordering).
            I = IAll(j,:);
            
            % Compute the corresponding marginal distribution (including
            % the target node itself).
            vv = margpdf(v,[i,I]);
            
            % Compute the error of all the functions.
            Err = sum(F.*(vv(1:kk)'*Ones1nf)) + ...
                sum((1-F).*(vv(kk+1:2*kk)'*Ones1nf));
            
            % Find the optimal (Bayes) k-variable predictor for the current
            % input variable combination.
            OptErr = min(Err);
            
            if OptErr<minerr
                % Select one of the optimal functions randomly.
                ind = find(Err==OptErr);
                ind = ind(randintex(1,1,length(ind)));
                Fhat(:,i) = F(:,ind);
                varF(:,i) = I';
                errF(i) = OptErr;
                minerr = OptErr;
            end % if OptErr<minerr
                        
        end % for i=starti:stopi
    end % for i=1:n
end % if isempty(F)
