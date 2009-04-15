function tf = ispercolating(I,c)

% tf = ispercolating(I,c) - check the percolation of a 2-D lattice
% This function checks that whether or not the given binary lattice (image)
% I is percolating. The image I is assumed to be binary valued where the
% pixels of interest are marked with ones and "background" pixels are set
% to zero. Current implementation uses so called "R1" rule when checking
% the percolation. That is, this function checks that whether or not there
% exists a cluster (connected component of ones) which spans the whole
% image in horizontal direction, i.e., a cluster touches both the left and
% right border of the image. The second parameter c defines the
% connectivity rule which must be either 4 or 8. In order to able to use
% this function one must have the Image Processing Toolbox.

% Functions used: bwlabel.m

% 01/10/2002 by Harri Lähdesmäki

w = size(I,2);              % Width of the lattice.

tf = 0;                     % Let's be conservative, the default value is False (0).

I = bwlabel(I,c);           % Find connected components, c specifies the connectivity.

% Check the span of connected components only to horizontal directions. That
% is, we assume "free boundaries to vertical directions". Check one connected 
% component at a time.
maxI = max(I(:));
for j=1:maxI
    [ind1,ind2] = find(I==j);
    if min(ind2)==1 & max(ind2)==w
        tf = 1;
        return;             % Exit whenever a percolating cluster is found.
    end
end
