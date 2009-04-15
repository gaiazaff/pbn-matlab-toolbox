function pbnDrawA(A)

% pbnDrawA(A) - visualize a state transition matrix
%
% This function *attempts* to draw the state-transition diagram for a
% BN/PBN with state transition matrix A. It uses an existing MATLAB package
% for drawing directed graphs. The placement of the states is determined by
% multidimensional scaling (MDS) where if state A has a high probability of
% going to state B, then, they are placed close to each other on the
% diagram. Of course, this doesn't always work well because the matrix A is
% not symmetric and the graph algorithm is quite primitive.

% Ilya Shmulevich; Aug. 28, 2001.
% Modified May 16, 2003 by HL. An updated Matlab package for drawing graphs
% is used. Package is located in a subfolder '\GraphLayout'.

path(path,'.\GraphLayout');             % add the drawgraph package to the path
nn = size(A,1);                         % number of states
n = log2(nn);                           % number of genes
%C = domds(A);                          % make sure domds is set up for 2 dimensions
[x,y,z] = cylinder(1,nn);
x = x(1,:);
y = y(1,:);
%x = C(:,1)';                           % get MDS coordinates in 2-D
%y = C(:,2)';                           % and store them into x and y
%x = x/range(x);                        % scale x and y to be between 0 and 1
%x = x - min(x);
%y = y/range(y);
%y = y - min(y);

for i = 1:nn,
   labels{i} = dec2bin(i-1,n);          % create labels
end

AA = zeros(size(A));
AA(find(A>0)) = 1;
%[dummyx, dummyy, h] = draw_layout(A>0,labels,zeros(1,nn),x,y);
[dummyx, dummyy, h] = draw_layout(AA,labels,zeros(1,nn),x,y);
axis auto

for i=1:length(h),
  col = rand(1,3);
  % patches
  set(h(i,2),'facecolor', col); drawnow;
  % text
  set(h(i,1),'color', 1-col); drawnow;
end;