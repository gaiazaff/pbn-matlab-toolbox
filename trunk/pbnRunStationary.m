function pmf = pbnRunStationary(x,F,varF,nf,nv,cij,p,nsteps)

% pmf = pbnRunStationary(x,F,varF,nf,nv,cij,p,nsteps)
% This function combines pbnRun and pbnStationary. The basic reason is that
% we may not want to store a huge "history" Y of the entire run. The
% parameters are the same as for pbnRun, but the output is the pmf vector,
% as in pbnStationary

% Functions used: pbnNextState.m

% Ilya Shmulevich; May 20, 2002
% Modified May 14, 2003 by HL.


n = length(nf);                                  % number of genes

if isstr(x) % if we want a random start
    x = rand(1,n)>0.5;
end

pmf = zeros(1,2^n);                             % initialize the pmf vector
pow = 2.^[n-1:-1:0]';                           % for conversion to decimal

for k=1:nsteps
    x = pbnNextState(x,F,varF,nf,nv,cij,p); % One step of the PBN.
    decx = x*pow + 1;                           % convert state to decimal
    pmf(decx) = pmf(decx) + 1;                  % update pmf at that state
end % for k=1:nsteps

pmf = pmf/nsteps;                               % then normalize
