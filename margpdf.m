function ff = margpdf(f,vars)

% ff = margpdf(f,vars) - joint marginal probability distribution
% Function returns a joint marginal probability distribution of f for the
% variables defined in vars. Note that the order of variables is
% significant. Vector f is assumed to represent a pdf of a binary-vector
% random variable. In detail, assume a binary vector random variable X =
% (X1,X2,...,Xn) having n elements. Then, f(1) = Pr{X = 00...00}, f(2) =
% Pr{X = 00...01}, ..., f(2^n) = Pr{X = 11...11}. vars is a vector having a
% subset of {1,...,n} as its elements. If vars is (Xi1,Xi2,...,Xik), then
% ff represents the joint marginal pdf for the corresponding variables (and
% in the same order).
% INPUT:
% f     - A row vector whose elements sum up to unity. The length of f must
%         be a power of two. 
% vars  - A vector of indices such that, when considered as a set, it is a
%         subset of {1,...,n}. Note that the ordering of the elements in
%         vector vars is relevant.
% OUTPUT:
% ff    - A vector representing the joint marginal pdf of f for the
%         variables defined in vars. 

% Harri Lähdesmäki 06/05/2003

n = round(log2(length(f))); % The number of variables in the random binary vector-variable.
nm = length(vars); % The number of variables in the joint marginal distribution.

ff = -ones(1,2^nm);

ind = zeros(1,2^n); % Initialize an index vector.
states = [0:(2^n-1)]; % All possible states of the orig. n variable random vector as a decimal number.
b = 2.^[(nm-1):-1:0];

% Element ind(i) is set to a decimal number based on the (binary) values of
% the marginal variables of the state i.
for i=1:nm
    ind = ind + bitget(states,n-vars(i)+1)*b(i);
end

%states = ff2n(n);
%ind = (states(:,vars)*b')';

% "Integrate" over the dummy variables.
for i=0:(2^nm-1)
    ff(i+1) = sum((ind==i).*f);
end
