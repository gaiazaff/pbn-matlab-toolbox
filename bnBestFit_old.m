function E = bnBestFit_old(D,IO,Ind,w,k,emax,pinfo)

% E = bnBestFit_old(D,IO,Ind,w,k,emax,pinfo) - Gene regulatory network inference
%
% bnBestFit performs a gene regulatory network inference under the Boolean
% network model for a set of "interesting" genes. The function outputs a
% "row-matrix" E containing the histograms of the error size (ref. Best-Fit
% Extension Paradigm) for all the interesting genes.
%
% Input variables:
% D     - Binary input matrix. D(i,:) corresponds to the (binary) expression
%         profile of the i:th gene. Correspondingly, D(:,j) represents the
%         gene activity profile for the j:th time index.
% IO    - Indices which define the input-output column pairs in the data
%         matrix D. The first and the second row define the input and output
%         columns in the matrix D. For example, IO(1,i) defines the i:th
%         input column vectors (the i:th gene activity profile), and IO(2,i)
%         defines the corresponding output activity profile. Each index in
%         IO must be between 1 and size(D,2). This parameter allows to
%         infer functions from separate time-series or even from the
%         steady-state data. That is, one can also set IO(1,i) == IO(2,i).
% Ind   - Index vector containing indices of the "interesting" genes. All
%         elements of Ind must be between 1 and size(D,1).
% w     - Weight vector containing positive weights for the measurements in
%         D. For now, the implementation only allows to define a single
%         weight for each column in D. In particular, the weight w(i)
%         defines the weight for the i:th input vector.
% k     - Maximum indegree.
% emax  - Maximum allowable error size.
% pinfo - An optional vector that defines the start and stop index of
%         the "main" loop (see below) and the first variable combination.
%         The start and stop indices are defined by pinfo(1) and pinfo(2).
%         The first variable combination is defined by the elements
%         pinfo(3:end) = pinfo(3:k+2)
%
% Output variables:
% E     - "Row-matrix" of the histograms of the error size (ref. Best-Fit
%         Extension Paradigm) for all the "interesting" genes.
%
% Functions used: BINARRAY2DEC, NEXTNCHOOSEK, DISTORTIONVECTORS

% 09.04.2002 by Harri Lähdesmäki
% Modified: 22.10.2002 by HL: The code was sped up. 29.10.2002 by HL:
% Name changed and one input parameter added. Now more than one time-series
% (i.e. concatenation of several time-series) can be used. Even a
% combination of time-series and steady-state data is allowed.

% Note that one possible way of speeding up the code is to form separate
% input and output data matrices. That would decrease the number of indexing
% in the main loop. This may, however, cause problems in the case of huge
% data sets and large number of "interesting" genes... Also, function
% DistortionVectors should be "checked".


% Two to the power of k (needed often).
kk = 2^k;

% The number of all genes.
n = size(D,1);
% The number of "useable" measurements.
m = size(IO,2);
%m = size(D,2)-1;
% The number of "interesting" genes.
ni = length(Ind);

%------------------------------
% Output values.
%Out = D(Ind,IO(2,:));

% Input vectors (examples).
%D = D(:,IO(1,:));
%------------------------------

% This matrix (C01) has the role of c^(0) and c^(1) for all interesting
% genes. Further, C01 = [c^(0),c^(1)];
C01 = zeros(ni,2*kk);
sC01 = size(C01);

% A temporary index vector. Not currently used. This is needed, however, if
% binary input strings are converted into decimal numbers by using a method
% based on bitset-function.
t1 = [1:ni];

% Powers of two (used in binarray2dec).
b = 2.^[k-1:-1:0]';

% Initialize the output to zero matrix.
E = zeros(ni,emax+1);

% The number of different variable combinations.
varnum = nchoosek(n,k);

% If pinfo exists, the start and stop indices of the "main" loop are taken
% from it.
% The first k variables, i.e., I = [1,2,...,k];.
I = [1:k];
starti = 1;
stopi = varnum;
if nargin>=7
    starti = pinfo(1);
    stopi = pinfo(2);
    I = pinfo(3:k+2);
end % if nargin>=7

disp('Process info:');
disp(['Start index: ',num2str(starti)]);
disp(['Stop index: ',num2str(stopi)]);
disp(['First variable comb.: ',num2str(I)]);

% Run through all variable combinations.
for i=starti:stopi
    
    % Initialize again.
    C01 = zeros(ni,2*kk);
    
    % Run through all measurements. This loop also takes into account possible
    % multiplisities in measurements, i.e., computes new weights for those
    % measurements that appear several times in T (and/or F).
    for j=1:m
        % The current input as a decimal number to be used to index the matrix
        % C01.
        %dn = binarr2dec(D(I,IO(1,j))',b) + 1;
        % Try to avoid unneccary function calls... does the same thing as the
        % previous line.
        dn = D(I,IO(1,j))'*b + 1;
        % Yet another way to convert the binary string into a decimal
        % number. Seems to be as fast as the previous line.
        %dn = sum(bitset(0,t1(logical(D(I,IO(1,j)))))) + 1;
        
        % Update C01 (c^(0) and c^(1)).
        % First update left half (C0) and then right half (C1)
        C01(logical(1-D(Ind,IO(2,j))),dn) = C01(logical(1-D(Ind,IO(2,j))),dn) + w(j);
        C01(logical(D(Ind,IO(2,j))),kk+dn) = C01(logical(D(Ind,IO(2,j))),kk+dn) + w(j);
    end % for j=1:m
    
    % Absolute difference between c^(0) and c^(1).
    C = abs(C01(:,1:kk)-C01(:,kk+1:end));
    
    % Find all proper functions using this variable combination. Perform the
    % search for all the interesting nodes.
    for j=1:ni
        % Compute the optimal function and the corresponding error size.
        t6 = [C01(j,1:kk);C01(j,kk+1:end)];
        t5 = min(t6);
        %[t5,fopt] = min(t6);
        % fopt = fopt - 1;
        % fopt = 1 - fopt;
        eopt = sum(t5);
        
        % Find proper distortion vector, if the optimal error size is less
        % than or equal to emax.
        if eopt<=emax
            % Sort |c^(0)-c^(1)|.
            c = sort(abs(C(j,:)));
            % The first distortion vector and set of proper distortion vectors.
            d = zeros(1,kk);
            DD = d;
            % Apply function DistortionVectors that will return a row-matrix of
            % proper (permuted) distortion vectors.
            DD = DistortionVectors(d,eopt,0,DD,c,emax);
            % Compute the error sizes for all proper (permuted) distortion
            % vectors.
            DD = eopt*ones(size(DD,1),1) + DD*c';
            % Update the distribution (histogram) of error sizes for the j:th
            % gene. For loop seems to be faster than the method based on
            % hist-function.
            %h = hist(DD,[0:emax]);
            %E(j,:) = E(j,:) + h;
            for jj=1:size(DD,1)
                E(j,DD(jj)+1) = E(j,DD(jj)+1) + 1;
            end % for jj=1:size(DD,1)
        end % if eopt<=emax
    end % for j=1:ni
    
    % The "next" variable combinations.
    I = nextnchoosek(I,n);
    
    % Display something.
    if mod(i,1000000)==0
        disp([num2str(i),'/',num2str(stopi)]);
        %save('resultstemp.mat','E','i','I');
    end % if mod(i,1000000)==0
    
end % for i=1:varnum

%save('resultstemp.mat','E','i','I');
