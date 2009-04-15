function tf = isAmuSmall(f,varargin)

% SYNTAX: tf = isAmuSmall(f,mu,[B])
% This function returns 1 or 0 (true or false) depending on whether the Boolean function 
% specified by the truth table f belongs on the class <a^mu>, where mu is an input 
% parameter. mu should also be at least 2.
% Ilya Shmulevich; 12/04/02

tf = isAmu(~flipud(f),varargin{:});
