function [p,c] = percolationcurve(Fall,sl)

% [p,c] = percolationcurve(Fall,[sl]) - percolation on a 2-D lattice
% This function creates the "percolation curve" for a specific class of
% Boolean functions defined by the input variable Fall, or for the class of
% all 4-variable Boolean functions if Fall is an empty matrix.
% INPUT:
% Fall  - contains all the specific functions used in the simulation (or is
%         empty). Each column of Fall is supposed to represent the truth
%         table of a 4-variable Boolean function.
% sl    - [Optional] defines the size of the square 2-D lattice.
% OUTPUT:
% The outputs p and c correspond to the x- and y-axes of the
% percolating curve. In the case of empty-matrix, p defines the bias of
% random boolean functions. Note that some of the parameter appearing in 
% the code below can also be adjusted even though they are not input 
% parameters for the function.

% Functions used: wiringsquare, bnRun, seqperiod, ispercolating

% 27/09/2002 by Harri Lähdesmäki, Modified: 01/10/2002 and May 14, 2003 by
% HL.

%-------------------------------
% User adjustable parameters.
%-------------------------------
maxiter = 100;              % Number of iterations in Monte-Carlo.
nsteps = 1500;              % Number of steps to run the network.
remsteps = 500;             % Number of "transient" steps to be removed.
p = [0:0.02:1];             % Selection probabilities for functions from Fall.
conrule = 4;                % Connectivity rule for the 2-D image.
%-------------------------------

bit = ~isempty(Fall);

if nargin<2
    sl = 50;                % Default size of the lattice (one dimension).
end

n = sl^2;                   % Number of genes.
nv = 4*ones(1,n);           % Number of variables per function.
lt = 2^4;                   % Default length of the truth table.

if bit==1                   % if Fall is non-empty.
    nf = size(Fall,2);      % Number of functions in the "set" Fall.
end

lp = length(p);             % Length of the p vector.

varF = wiringsquare(sl);    % Square wiring for 2-D lattice, fixed.
F = zeros(lt,n);            % Allocate memory for the random truth tables (to be generated later).
c = zeros(1,lp);            % Allocate memory for the output.

for k=1:lp                  % Do the Monte-Carlo for all p's.    
    for i=1:maxiter         % For each element of p, generate maxiter number of random networks.
        
        if bit
            % Select (randomly) the nodes that are going to have a randomly chosen function from Fall.
            ind = rand(1,n)>(1-p(k));
            % Number of functions selected from Fall.
            sumoffind = sum(ind);
            F(:,find(ind==1)) = Fall(:,unidrnd(nf,1,sumoffind));
            % Other functions are selected uniformly randomly.
            F(:,find(ind==0)) = rand(lt,n-sumoffind)>0.5;
        else
            % Select uniformly randomly the genes that are biased up and down.
            ind = rand(1,n)>0.5;
            F = rand(lt,n);
            F(:,find(ind==0)) = F(:,find(ind==0))>(1-p(k));     % Bias up.
            F(:,find(ind==1)) = F(:,find(ind==1))>p(k);         % Bias down.
        end
        
        % Run the Boolean network.
        Y = bnRun('rand',F,varF,nv,remsteps+1);
        Y = bnRun(Y(end,:),F,varF,nv,nsteps-remsteps-1);
        % Remove the first "removesteps" steps. (Note that the first step is the 'random' point.)
        %Y = Y(remsteps+2:end,:);
                
        Y = seqperiod(Y);       % Minimum-length repeating sequence for each gene.
        Y = reshape(Y,sl,sl)';  % Reshape the period lengths into a (2-D) matrix.
        Y = (Y==1);             % Find the elements that are fized, i.e., have period length 1.
        
        % Check the percolation and add the counter in the affirmative case.
        c(k) = c(k) + ispercolating(Y,conrule);
        
    end
    disp([num2str(k),'/',num2str(lp)]);
end

c = c/maxiter;              % Normalize the counts.
