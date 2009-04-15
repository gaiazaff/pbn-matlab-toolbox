function y = pbnRun(x,F,varF,nf,nv,cij,p,nsteps)

% y = pbnRun(x,F,varF,nf,nv,cij,p,nsteps) - run a PBN
% This function runs a PBN starting from the initial starting state x for
% nsteps steps. It essentially uses pbnNextState in an iterative fashion.
% The output is a matrix Y (nsteps+1-by-length(x)) containing the history
% of the process. If x is not a vector, but the word 'rand', then a random
% starting state is used

% Functions used: pbnNextState.m

% Ilya Shmulevich; Aug. 14, 2001
% Modified May 12, 2003 by HL.

n = length(nf); % number of genes

if isstr(x) % if we want a random start
    x = rand(1,n)>0.5;
end

y = zeros(nsteps+1,n); % initialize y
y(1,:) = x; % first one is initial state
for step=2:nsteps+1,
   x = pbnNextState(x,F,varF,nf,nv,cij,p); %otherwise use logical functions
   y(step,:) = x; % update history
end
