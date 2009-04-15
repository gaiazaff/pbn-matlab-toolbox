function y = randintex(n,low,high)

% y = randintex(n,low,high) - n random non-repeated integers between low and high
% randintex(n,low,high) returns n uniformly random integers between low and
% high (limits included) taking into account and preventing repetition.
% Restrictions: low<=high and 0<=n<=high-low+1.

% Functions used:

% 15.08.2001 by Harri Lähdesmäki ©
% Modified: 30.09.2002 by HL: Complete new implementation.

% Check the restrictions.
if low>high
   error('low must satisfy low<=high.');
end % if low>high
if n>(high-low+1) | n<0
   error('n must satisfy 0<=n<=high-low+1.');
end % if n>(high-low+1) | n<0

% Number of integers between (and including) low and high.
N = high-low+1;

y = zeros(1,n);             % Allocate memory for the output.

integers = [low:high];

for i=1:n
    ind = unidrnd(N+1-i);   % Select an integer from [1,..,N+1-i] uniformly randomly.
    y(i) = integers(ind);   % The element ind of Integers is one of the output values.
    integers(ind) = [];     % Remove the selected integer so that it can not be selected anymore.
end

