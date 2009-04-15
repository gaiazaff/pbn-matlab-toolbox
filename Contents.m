function Contents

% BN/PBN toolbox.
% May 16, 2003.
%
% The way of representing a BN/PBN in this toolbox is defined e.g. in
% pbnRnd.m
%
% RUN NETWORKS.
%
% bnNextState - one step of a Boolean Network
% bnRun - run a Boolean Network
% pbnNextState - one step of a PBN
% pbnRun - run a PBN
%
%
% NETWORK STATISTICS.
%
% bnDerridaCurve - Derrida curve
% bnStats - statistics of a Boolean Network
% ispercolating - check the percolation of a 2-D lattice
% pbnInfluence - Influence matrix of a PBN
% percolationcurve - percolation on a 2-D lattice
% wiringsquare - wiring for a square-lattice Boolean network
%
%
% STATE TRANSITIONS & DISTRIBUTIONS.
%
% bnA - state transition matrix a Boolean network
% bnAsparse - sparse state transition matrix a Boolean network
% bnAttractor - attractors and distances to the attractors
% bnAttractorMarkov - attractor proximity matrix
% pbnA - state transition matrix of a PBN
% pbnAij - transition prob. between two states
% pbnRunStationary - combines pbnRun and pbnStationary
% pbnStationary - empirical stationary distribution
%
%
% INFERENCE.
%
% bnBestFit - Best-Fit inference
% bnBestFit_old - Gene regulatory network inference
% bnCrossVal - Cross-validation for Boolean Network inference
% pbnBayesPred - Bayes predictors for a PBN
%
%
% RANDOM QUANTITIES.
%
% margpdf - joint marginal probability distribution
% pbnRnd - random BN/PBN
% pmfrnd - sample a discrete distribution
% randintex - n random non-repeated integers between low and high
% rndA2
% rndcanalnew
% rndperm - random permutation of a set of integers
%
%
% VISUALIZATION & PRINTING.
%
% pbnDrawA - visualize a state transition matrix
% printPost
% SquareRun
%
%
% INTERVENTION.
%
% pbnFindFunc - alter the functions of a PBN
% pbnStateVisit - best gene for intervention
%
%
% MEMBERSHIP TESTING & BOOLEAN FUNCTIONS.
%
% buildB
% isA5
% isAinf
% isAinfSmall
% isAmu
% isAmuSmall
% iscanal
% isselfdual
% permutevars - Permute variables in a Boolean function
%
%
% MISC.
%
% makeK - matrix K
