function tf = isselfdual(f)

% SYNTAX: tf = isselfdual(f)
% This function takes a column vector representing a truth table of a Boolean function
% and determines whether or not this function is self-dual.
% Ilya Shmulevich; 08/23/02
% Modified: 24/02/2003 by IS.

tf = 0;                         % default is FALSE
n = length(f)/2;                % half of the length of the truth table

%tf = all([eye(n) rot90(eye(n))]*f);
%tf = (sum([eye(n) rot90(eye(n))]*f) == n);
tf = ([eye(n) rot90(eye(n))]*f) == ones(n,1);
