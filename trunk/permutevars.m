function ff = permutevars(f,vars,varsnew)

% ff = permutevars(f,vars,varsnew) - Permute variables in a Boolean function
% This function takes as an input a Boolean function f, indices of its
% current input variables (in vars), and the possibly a different
% permutation of the variables in varsnew. The function returns the truth
% table of the same function corresponding to the new variable ordering.
% INPUT:
% f         - A truth table of a Boolean function.
% vars      - Indices of input variables of the function f.
% varsnew   - Possibly a different ordering of the input variables.
% OUTPUT:
% ff        - Truth table of the function f corresponding to the new
%             variable ordering.

% May 15, 2003 by Harri Lähdesmäki.

k = length(vars); % The number of variables.
b = 2.^[k-1:-1:0]';
x = ff2n(k); % All possible input vectors.

ind = zeros(1,k);
for i=1:k
    ind(i) = find(vars==varsnew(i)); % Position of varsnew(i) in vars.
end % for i=1:k

x = x(:,ind); % Permute the columns of x according to ind.
x = x*b + 1; % Convert the binary numbers (rows of x) into a decimal number.
ff = zeros(2^k,1);
ff(x) = f; % Permute f according to x.
