function [Fhat,OptE,Et] = bnBestFit(X,Y,w,k,F,Xt,Yt,wt)

% [Fhat,OptE,Et] = bnBestFit(X,Y,w,k,F,Xt,Yt,wt) - Best-Fit inference
%
% bnBestFit performs a gene regulatory network inference under the Boolean
% network model for a set of genes. The function returns the Best-Fit
% function and the corresponding error-size for all input variable
% combinations (in X) and for all the genes in Y. That is, rows in X
% correspond to predictor genes, and rows in Y correspond to target genes.
% Currently, the i:th gene at the j:th sample (time point) Y(i,j) is
% predicted based on the predictor gene values at the same sample X(:,j).
% Note that if unity weights (defined in w) are used for all the samples 
% the found functions are equal to the ones which minimize the
% resubstitution error on the sample data. In such a case, the
% corresponding error-sizes in OptE are the (non-normalized) minimum
% resubstitution errors, i.e., the number of error/misclassifications. In
% case of tie, a random selected function with minimum error-size is
% returned. If additional test data sets (Xt, Yt, and wt) are provided the
% error-size of the found Best-Fit functions on the test data is compute as
% well. (This can be useful in the case of cross-validation experiments.)
%
% INPUT:
% X     - Binary input matrix. X(i,:) corresponds to the (binary)
%         expression profile of the i:th predictor gene. Correspondingly,
%         X(:,j) represents the gene activity profile of the predictor
%         genes for the j:th sample.
% Y     - Binary output matrix. Y(i,:) corresponds to the (binary)
%         expression profile of the i:th target gene. Correspondingly,
%         Y(:,j) represents the gene activity profile of the target genes
%         for the j:th sample. Y(i,j) is the value of the i:th gene for the
%         j:th sample, i.e., the bit that is to be predicted based on
%         X(:,j). (Holds for all i.)
% w     - Weight vector containing positive weights for the measurements in
%         X and Y. For now, the implementation only allows to define a
%         single weight for each column in X (and Y). In particular, the
%         weight w(i) defines the weight for the i:th input vector (and the
%         corresponding output).
% k     - Maximum indegree of the predictor functions.
% F     - The set of Boolean predictors to be used in the inference. F is a
%         (2^k)-by-nf binary matrix, where k is the number of variables in
%         each function and nf is the number of functions. If F is an empty
%         matrix, then the function class is considered to contain all
%         k-variable Boolean functions. Let f = F(:,j) be the j:th column
%         of F. Then, f(0) defines the output value for input vector
%         00...00, f(1) for 00...01, f(2) for 00...10, ..., and f(2^k) for
%         11...11. Input vectors are interpreted such that the k:th (left
%         most) bit defines the value of the first input variable, the
%         k-1:th bit defines the value of the second input variable, ...,
%         and the right most (first bit) defines the value of the last
%         input variable.
% Xt    - [Optional] Input data for a separate test data. Format is the
%         same as for the matrix X (see above).
% Yt    - [Optional] Output data for a separate test data. Format is the
%         same as for the matrix Y (see above).
% wt    - [Optional] Weights for a separate test data. Format is the same
%         as for the matrix w (see above).
%
% OUTPUT:
% Fhat  - A 3-D matrix of the Best-Fit functions for each input variable
%         combinations and for all nodes. Fhat has size
%         (2^k)-by-nchoosek(n,k)-by-ni, where n is the size of the first
%         dimension of matrix X (i.e., the number of predictor genes) and
%         ni is the number of target nodes. Fhat(:,:,i) defines the
%         Best-Fit functions for the i:th node. Each column in Fhat(:,:,i)
%         is interpreted as the columns in F (see above).
% OptE  - The error-size of the Best-Fit function for all input
%         combinations and for all the target nodes. OptE has size
%         nchoosek(n,k)-by-ni. Thus, OptE(i,j) is the error-size of the
%         Best-Fit function for i:th input variable combination and for the
%         j:th target node.
% Et    - This variable is returned only if Xt, Yt and wt are present in
%         the input. The error-size of the Best-Fit function for all input
%         combinations and for all the target nodes on a separate test
%         data.

% 03.04.2003 by Harri Lähdesmäki, modified from bnBestFit.
% Modified May 14, 2003 by HL.


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define and initialize some variables.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

[n,m] = size(X); % The number of predictor genes and the number of measurements.
ni = size(Y,1); % The number of (target) genes.

b = 2.^[k-1:-1:0]'; % Powers of two (used in binary-to-decimal convertions).

W = ones(ni,1)*w; % Weights in a matrix form (assume w is a row vector).
%W = repmat(w,ni,1); % The same as above.

kk = 2^k; % Two to the power of k (needed often).

% A temporary index vector. Not currently used. This is needed, however, if
% binary input strings are converted into decimal numbers by using a method
% based on bitset-function.
%t1 = [1:ni];

combnum = nchoosek(n,k); % The number of different variable combinations.

% Generate all variable combinations in advance. This will work only for
% moderately small data sets. If one wants to use larger data sets, then
% the input variable combinations can be generated using the function
% nextnchoosek.m (e.g. given an input variable combination, I =
% nextnchoosek(I,n); generates the next variable combination in
% lexicographial order.
if combnum>20000 % Limit the number of possible combinations.
    error('Too many variable combinations. Modify the code a little bit...')
end % if combnum>20000
IAll = nchoosek([1:n],k);

% Modify the variables below if only a subset of all combinations are to be
% checked (e.g. in the case of parallelizing the code...)
starti = 1;
stopi = combnum;

% Initialize the output matrix/matrices.
Fhat = zeros(kk,combnum,ni);
OptE = zeros(combnum,ni);

TestBit = 0;
if nargin==8 % If the separate test data sets are present in the input.
    Et = zeros(combnum,ni);
    Wt = ones(ni,1)*wt; % Wt = repmat(w,ni,1);, weights in vector form.
    TestBit = 1;
end % if nargout > 1



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The main loop separately for unconstrained and constrained case.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% If F is an empty matrix, then the function class is considered to contain
% all k-variables Boolean functions (unconstrained case).
if isempty(F)
    
    % Two times two to the power of k (needed often).
    kkk = 2*kk;
    
    % This matrix (C01) has the role of c^(0) and c^(1) for all interesting
    % genes. Further, C01 = [c^(0),c^(1)];
    C01 = zeros(ni,kkk);
    %sC01 = size(C01);
    
    % Indices of the undefined inputs (not occuring in the training data).
    %Undef = ones(kk,1);
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Run through all variable combinations.
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for i=starti:stopi
        
        % The current variable combinations (in lexicographical ordering).
        I = IAll(i,:);
        
        % Initialize again.
        C01 = zeros(ni,kkk);
        %Undef = ones(kk,1);
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Run through all measurements.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % This loop also takes into account possible multiplisities in
        % measurements, i.e., computes new weights for those measurements
        % that appear several times in T (and/or F).
        for j=1:m
            
            % The current input as a decimal number to be used to index the
            % matrix C01.
            %dn = binarr2dec(D(I,IO(1,j))',b) + 1;
            
            % Try to avoid unneccary function calls... does the same thing
            % as the previous command.
            dn = X(I,j)'*b + 1;
            
            % Yet another way to convert the binary string into a decimal
            % number. Seems to be approx. as fast as the previous command.
            %dn = sum(bitset(0,t1(logical(D(I,IO(1,j)))))) + 1;
            
            % Update C01 (c^(0) and c^(1)). First update left half (C0) and
            % then right half (C1)
            C01(logical(1-Y(:,j)),dn) = C01(logical(1-Y(:,j)),dn) + w(j);
            C01(logical(Y(:,j)),kk+dn) = C01(logical(Y(:,j)),kk+dn) + w(j);
            
            % Update Undef.
            %Undef(dn) = 0;
        end % for j=1:m
        
        % Find the Best-Fit function for all the nodes.
        [OptErr,OptF] = min(cat(3,C01(:,1:kk),C01(:,kk+1:end)),[],3);
        %OptF = OptF' - 1; % Change the scale from {1,2} to {0,1}.
        %Fhat(:,i,:) = 1 - OptF; % Change the bits.
        %Fhat(:,i,:) = 2 - OptF'; % The same as above.
        OptF = 2 - OptF';
        
        % Below, all undefined bits in the truth tables are flipped
        % uniformly randomly.
        %OptF(logical(Undef),:) = ((rand(sum(Undef),ni))>0.5);
        
        % All output bits having tie are set uniformly randomly. This also
        % takes care of the undefined bits due to the initialization of
        % matrix C01.
        Ties = (C01(:,1:kk)==C01(:,kk+1:end))';
        OptF(Ties) = (rand(1,sum(Ties(:))))>0.5;
        
        % Store the Best-Fit functions.
        Fhat(:,i,:) = OptF;
        % Store the corresponding (weighted) error-size.
        OptE(i,:) = sum(OptErr,2)';
        
        
%        % Modify the code below if all functions with limited error-size
%        % should be found.
%         %================================================================
%         % Absolute difference between c^(0) and c^(1).
%         C = abs(C01(:,1:kk)-C01(:,kk+1:end));
%         
%         % Find all proper functions using this variable combination. Perform the
%         % search for all the interesting nodes.
%         for j=1:ni
%             % Compute the optimal function and the corresponding error size.
%             t6 = [C01(j,1:kk);C01(j,kk+1:end)];
%             t5 = min(t6);
%             %[t5,fopt] = min(t6);
%             % fopt = fopt - 1;
%             % fopt = 1 - fopt;
%             eopt = sum(t5);
%             
%             % Find proper distortion vector, if the optimal error size is less
%             % than or equal to emax.
%             if eopt<=emax
%                 % Sort |c^(0)-c^(1)|.
%                 c = sort(abs(C(j,:)));
%                 % The first distortion vector and set of proper distortion vectors.
%                 d = zeros(1,kk);
%                 DD = d;
%                 % Apply function DistortionVectors that will return a row-matrix of
%                 % proper (permuted) distortion vectors.
%                 DD = DistortionVectors(d,eopt,0,DD,c,emax);
%                 % Compute the error sizes for all proper (permuted) distortion
%                 % vectors.
%                 DD = eopt*ones(size(DD,1),1) + DD*c';
%                 % Update the distribution (histogram) of error sizes for the j:th
%                 % gene. For loop seems to be faster than the method based on
%                 % hist-function.
%                 %h = hist(DD,[0:emax]);
%                 %E(j,:) = E(j,:) + h;
%                 for jj=1:size(DD,1)
%                     E(j,DD(jj)+1) = E(j,DD(jj)+1) + 1;
%                 end % for jj=1:size(DD,1)
%             end % if eopt<=emax
%         end % for j=1:ni
%         %================================================================
        
        
        if TestBit % If the test data is provided
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            % Apply the Best-Fit functions to the test data.
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            % All the inputs on the test data as decimal number.
            dn = Xt(I,:)'*b + 1;
            
            % Output values of the current functions for all the inputs.
            %yy = squeeze(Fhat(dn,j,:))';
            yy = OptF(dn,:)';
        
            % Compare the computed output to the desired one. Note that this
            % also takes into account possible multiplisities in measurements,
            % i.e., computes new weights for those measurements that appear
            % several times in T (and/or F).
            Et(i,:) = sum((yy~=Yt).*Wt,2)';
            %Et(i,:) = sum((yy~=Yt),2)';
        end % if TestBit
        
        
        % The next variable combinations (in lexicographical ordering)
        %I = nextnchoosek(I,n);
        
        % Display something.
        %if mod(i,1000)==0
        %    disp([num2str(i),'/',num2str(stopi)]);
        %    %save(file,'Fhat','i','I');
        %end % if mod(i,1000)==0
        
    end % for i=starti:stopi
    
else
    % If F is non-empty matrix, then the function class is defined in F,
    % i.e., the inference below is constrained by F.
    
    nf = size(F,2); % The number of functions in F.
    maxerror = sum(w) + 1; % The largest possible error plus one.
    ones_ni_1 = ones(ni,1);
    ones_1_nf = ones(1,nf);
    ones_ni_1_minus = -ones(ni,1);
    
    % Weights in a matrix form (assume w is a row vector).
    Wc = w'*ones(1,nf);
    %Wc = repmat(w',1,nf);
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Run through all variable combinations.
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for i=starti:stopi
        
        % In order to prevent a possible systematic bias caused by a
        % specific ordering of the functions in F in the case of tie, the
        % columns are randomly reorded for each variable combinations.
        F = F(:,randintex(nf,1,nf));
        
        % The current variable combinations (in lexicographical ordering).
        I = IAll(i,:);
        
        % Keep track of the minimum error-size and the corresponding
        % function (represented as an integer, which is a column number of
        % F).
        %bestesize = maxerror*ones_ni_1;
        %bestf = ones_ni_1_minus;
        
        % All the inputs as decimal number.
        dn = X(I,:)'*b + 1;
        
        % Output values of all the functions for all the inputs.
        yy = F(dn,:);
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Run through all the target genes.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Compare output values of all the functions to the desided output
        % of a single gene.
        for j=1:ni
            
            % Compute the error-size of all the function for the j:th node.
            E = (yy~=(Y(j,:)'*ones_1_nf)).*Wc;
            esize = sum(E);
            % Find the smallest error-size.
            [emin,indmin] = min(esize);
            % Store the the best candidate function.
            Fhat(:,i,j) = F(:,indmin);
            OptE(i,j) = emin;
            
        end % for j=1:ni
        
        
        
        if TestBit % If the test data is provided
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            % Apply the Best-Fit functions to the test data.
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            % All the inputs on the test data as decimal number.
            dn = Xt(I,:)'*b + 1;
            
            % Output values of the current functions for all the inputs.
            yy = squeeze(Fhat(dn,i,:))';
            
            % Compare the computed output to the desired one. Note that this
            % also takes into account possible multiplisities in measurements,
            % i.e., computes new weights for those measurements that appear
            % several times in T (and/or F).
            Et(i,:) = sum((yy~=Yt).*Wt,2)';
        end % if TestBit

        
%        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%        % Run through all the functions in F.
%        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%         % Run through all functions (exhaustive search).
%         for r=1:nf
%            
%             % Output values of the current function for all the inputs.
%             yy = F(dn,r)'; % row vector
%             %yy = repmat(yy,ni,1);
%             yy = ones_ni_1*yy; % Does the same thing as the row above.
%             
%             % Compare the computed output to the desired one. Note that
%             % this also takes into account possible multiplisities in
%             % measurements, i.e., computes new weights for those
%             % measurements that appear several times in T (and/or F).
%             E = (yy~=Y).*W;
%             %E = (mod(yy+Y,2)).*W;
%             esize = sum(E,2); % column vector
%             
%             % Find the errorsizes which are smaller than the previously
%             % found ones and update the best error-sizes as well as the
%             % indices of the best functions.
%             ind = (esize<bestesize);
%             bestesize(ind) = esize(ind);
%             bestf(ind) = r;
%             
%         end % for r=1:nf
%         
%         % Update the output matrix.
%         Fhat(:,i,:) = F(:,bestf);
%         if OutputBit
%             OptE(i,:) = bestesize';
%         end % if OutputBit
        
        % The next variable combinations (in lexicographical ordering)
        %I = nextnchoosek(I,n);
        
        % Display something.
        %if mod(i,100)==0
        %    disp([num2str(i),'/',num2str(stopi)]);            
        %    %save('resultstemp.mat','E','i','I');
        %end % if mod(i,1000000)==0
                
    end % for i=starti:stopi
end % if isempty(F)

%save(file,'E','i','I');
