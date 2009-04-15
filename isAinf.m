function tf = isAinf(f,A)

% SYNTAX: tf = isAinf(f,[A])
% This function takes a column vector representing a truth table of a Boolean function
% and determines whether or not this function belongs on the class <A^Inf>. The second 
% optional parameter (A) is the matrix that can be passed to speed things up, so it 
% doesn't have to be generated every time.
% Ilya Shmulevich; 08/22/02 and Harri Lahdesmaki 19/09/02

% Implementation is based on the function iscanal.m.

tf = 0;                         % default is FALSE
n = log2(length(f));            % number of variables

if nargin == 1                  % if A was not passed, let's create it
    A = ~ff2n(n)';              % A equals to the lower part of the transform matrix for testing forcing functions
end

c = A*(1-f);                    % take the transform. Note that only the ~f is transformed (see iscanal.m)
if any(c==2^(n-1))              % if c contains a value equal to one-half of the length of the truth table
    tf = 1;                     % then it's <A^Inf>
end