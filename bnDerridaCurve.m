function [dt,dt1] = bnDerridaCurve(F,varF,nv,nsteps)

% [dt,dt1] = bnDerridaCurve(F,varF,nv,nsteps) - Derrida curve
% This function creates the Derrida curves for a Boolean network. The
% network is specified by F, varF and nv (see e.g., bnNextState) nsteps
% specifies how many different random initial state pairs should be chosen.
% The outputs dt and dt1 correspond to the x- and y-axes of the Derrida
% curve.

% Functions used: randintex.m bnNextState.m

% Ilya Shmulevich; 09/30/02
% Modified May 14, 2003 by HL.

ngenes = size(F,2);
dt1 = zeros(1,ngenes);

for hd = 1:ngenes,                              % for each possible Hamming distance between starting states
    for i = 1:nsteps,
        x = unifrnd(0,1,1,ngenes)>0.5;          % choose initial state randomly
        y = x;                                  % let the other state equal x
        ind = randintex(hd,1,ngenes);           % Indices of the flipped nodes.
        y(ind) = 1 - y(ind);                    % Flip the selected bits.
        xx = bnNextState(x,F,varF,nv);          % Successor states.
        yy = bnNextState(y,F,varF,nv);
        
        dt1(hd) = dt1(hd) + sum(xx ~= yy);
    end
end

dt1 = (dt1 / nsteps) / ngenes;                  % average dt1 by dividing by number of steps and normalizing Hamming distance
dt = [1:ngenes]/ngenes;                         % create x-axis of Derrida curve