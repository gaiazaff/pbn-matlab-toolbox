function printPost(f)

% SYNTAX: printPost(f)
% This function does not output anything, but simply prints out all
% of the classes (A^mu, a^mu) to which the function f
% belongs. It uses isAmu and related functions.
% Note: if f is a member of A^mu and mu=n, where n is the number of
% variables, then f is also a member of A^Inf. Similarly for a^mu.
% Ilya Shmulevich; 02/04/02

f = f(:);
n = log2(length(f));

for mu = 2:n,
    if isAmu(f,mu)
        fprintf(1,'A^%d ',mu);
    end
    if isAmuSmall(f,mu)
        fprintf(1,'a^%d ',mu);
    end
end
fprintf(1,'\n');