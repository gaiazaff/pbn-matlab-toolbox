function varF = wiringsquare(n)

% varF = wiringsquare(n) - wiring for a square-lattice Boolean network
% This function outputs the wiring for a square-lattice Boolean network. The
% lattice is of size n x n and so, there are n^2 elements. varF is a 4xn^2
% matrix, in the same format as that used for the pbnA or bnA functions.
% The order of the variables (neighbors) is clockwise [left up right down]'
% This wiring assumes the 'wrap-around' boundary extension. So, the lattice
% is really a doughnut.

% Ilya Shmulevich; 08/28/02

varF = zeros(4,n^2);

for i = 1:n,
    for j = 1:n,
        if i == 1,
            up = n;
            down = i+1;
        elseif i == n
            down = 1;
            up = i-1;
        else
            up = i-1; down = i+1;
        end
        if j == 1
            left = n;
            right = j+1;
        elseif j == n
            right = 1;
            left = j-1;
        else
            left = j-1; right = j+1;
        end
        varF(:,sub2ind([n n],i,j)) = [sub2ind([n n],i,left) sub2ind([n n],up,j) sub2ind([n n],i,right) sub2ind([n n],down,j)]';
    end
end