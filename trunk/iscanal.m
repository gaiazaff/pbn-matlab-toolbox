function tf = iscanal(f,A)

% SYNTAX: tf = iscanal(f,[A])
% This function takes a column vector representing a truth table of a Boolean function
% and determines whether or not this function is canalizing.
% The second optional parameter (A) is the matrix that can be passed to speed things up, so it doesn't have to be generated every time.
% Ilya Shmulevich; 08/22/02

tf = 0;                         % default is FALSE
n = log2(length(f));            % number of variables

if nargin == 1                  % if A was not passed, let's create it
    T = ff2n(n)';
    A = [T ; ~T];
end

c = A*[f ~f];                   % take the transform
if any(any((c==2^(n-1))))       % if it contains a value equal to one-half of the length of the truth table
    tf = 1;                     % then it's canalizing
end