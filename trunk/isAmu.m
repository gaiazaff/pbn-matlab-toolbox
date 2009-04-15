function tf = isAmu(f,mu,B)

% SYNTAX: tf = isAmu(f,mu,[B])
% This function returns 1 or 0 (true or false) depending on whether the Boolean function 
% specified by the truth table f belongs on the class <A^mu>, where mu is an input 
% parameter. mu should also be at least 2.
% Ilya Shmulevich; 08/29/02
%
% Modified: 19/09/02 by HL
% Function was revised and renamed. A bug was corrected. Now it also works for mu = Inf. 
% However, it is recommended to use function isAinf for testing <A^Inf>. Matrix B is 
% assumed to be build such that it corresponds to the input parameter mu.

if f(1)
    tf = 0;
else
    n = log2(length(f));            % number of variables
    m = sum(f);                     % number of true vectors
    k = min([mu,m,n]);              % we only need to test <A^k>
    
    if nargin == 2 | k < mu         % if B was not passed, let's create it
        B = buildB(n,k);            % NOTE, build Bnk
    end
    
    c = B*f;
    tf = ~any(c == k);
end