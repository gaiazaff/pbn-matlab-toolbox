function I = nextnchoosek(I,N)

% I = nextnchoosek(I,N) - Next nchoosek variable combination
%
% nextnchoosek -function returns the "next" nchoosek combination of integers
% between 1 and N. The "next" integer combination is taken in lexicographical
% order. Note that this functions does not check that whether or not the
% input variable combination I is proper/valid.
%
% Input variables:
% I     - Previous integer combination.
% N     - The largest integer that can appear in I.
%
% Output variables:
% I     - The next variable combination or empty matrix if
%         I = [N-lenght(I)+1,...,N-2,N-1,N] is in the input.
%
% Functions used:

% 09.04.2002 by Harri Lähdesmäki
% Modified:

% Upper bounds for each element in I.
N = [N-length(I)+1:N];

% Start to check I from right to left that whether or not its elements can be
% increase by one.
for j=length(I):-1:1
    % Take action if j:th element can be increased by one.
    if I(j)<N(j)
        I(j) = I(j) + 1;
        p = I(j) + 1;
        % Change the elements that are right from the j:th position and quit.
        for k=(j+1):length(I)
            I(k) = p;
            p = p + 1;
        end
        return;
    end % if I(j)<N(j)
end % for j=length(I):-1:1

% Set the output to empty matrix if the input arguments were not proper.
I = [];
