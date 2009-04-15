function y = bnRun(x,F,varF,nv,nsteps)

% y = bnRun(x,F,varF,nv,nsteps) - run a Boolean Network
% This function runs a BN starting from the initial starting state x (a
% binary row vector) for nsteps steps. It essentially uses bnNextState in
% an iterative fashion. The output is a matrix Y (nsteps+1-by-length(x))
% containing the history of the process. If x is not a vector, but the word
% 'rand', then a random starting state is used

% Functions used: bnNextState.m

% Ilya Shmulevich; Aug. 14, 2001
% Modified May 12, 2003 by HL.

n = size(varF,2); % number of genes

if isstr(x) % if we want a random start
    x = rand(1,n)>0.5;
end

y = zeros(nsteps+1,n); % initialize Y
y(1,:) = x; % first one is initial state
for step = 2:nsteps+1,
   x = bnNextState(x,F,varF,nv);
   y(step,:) = x; % update history
end
