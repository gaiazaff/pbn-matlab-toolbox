function M=SquareRun(Y)

% SYNTAX: SquareRun(Y)
% This function, which has no outputs, takes a history Y and displays it as a square image.
% It is assumed that Y is of an appropriate length, so that it can be reshaped into a square image.
% Ilya; Aug. 28, 2002

n = sqrt(size(Y,2));
imshow(reshape(Y(1,:),n,n));
iptsetpref('ImshowTruesize','manual')
truesize([500 500])

for i = 2:size(Y,1),
    imshow(reshape(Y(i,:),n,n))
    drawnow
    %M(i-1) = getframe;
end
%movie(M)