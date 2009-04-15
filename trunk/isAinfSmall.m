function tf = isAinfSmall(f,varargin)

% SYNTAX: tf = isAinf(f,[A])
% This function takes a column vector representing a truth table of a Boolean function
% and determines whether or not this function belongs on the class <a^Inf>. The second 
% optional parameter (A) is the matrix that can be passed to speed things up, so it 
% doesn't have to be generated every time. A can be generated as: A = ~ff2n(n)'; 
% where n is the number of variables of the Boolean function.
% Ilya Shmulevich; 02/04/02

tf = isAinf(~flipud(f),varargin{:});