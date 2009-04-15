function Ehat = bnCrossVal(X,Y,w,k,F,cvk,nr)

% Ehat = bnCrossVal(X,Y,w,k,F,cvk,nr) - Cross-validation for Boolean Network inference
%
% Function estimates the (weighted) error of all predictor gene (rowes in
% X) combinations for all the target genes (rows in Y) using
% cross-validation. Currently, the predictor function inference itself is
% done using the bnBestFit.m function (other functions can be used as
% well). Note that if unity weights (defined in w) are used for all the
% samples , then  the estimated C-V error is equal to the standard error
% estimate.
%
% INPUT:
% X     - Binary input matrix. X(i,:) corresponds to the (binary)
%         expression profile of the i:th gene. Correspondingly, X(:,j)
%         represents the gene activity profile for the j:th sample.
% Y     - Binary output matrix. Y(i,j) is the output value of the i:th gene
%         for the j:th sample, i.e., the bit that is to be predicted based
%         on X(:,j). (Holds for all i.)
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
%         11...11. Input vectors are interpreted such that the first bit
%         defines the value of the first input node, the second bit defines
%         the value of the second input node, and so on.
% cvk   - The number of folds in the cross-validation. Limitation:
%         2<=cvk<=n.
% nr    - The number of times a single cross-validation procedure is
%         repeated.
%
% OUTPUT:
% Ehat  - Estimated error for all predictor gene combinations and for all
%         target genes. Ehat has size nchoosek(n,k)-by-ni, where n is the
%         number of predictor genes, k is the number of variables in the
%         Boolean functions, and ni is the number of target genes.

% Functions used: bnBestFit

% 03.04.2003 by Harri Lähdesmäki, modified from bnBestFit.
% Modified May 14, 2003 by HL.


% The number of genes and samples.
[n,m] = size(X);

% The number of target nodes.
ni = size(Y,1);

% Powers of two (used in binarray2dec).
b = 2.^[k-1:-1:0]';


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The number of samples in each split.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The number of samples in each split of the data.
ns = round(m/cvk)*ones(1,cvk);
ind = randintex(abs(m-sum(ns)),1,cvk);
ns(ind) = ns(ind) + sign(m-sum(ns));

% All input combinations.
combnum = nchoosek(n,k);
IAll = nchoosek([1:n],k);
Ehat = zeros(combnum,ni);
indAll = [1:m];


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The main loop.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Repeat the standard cross-validation nr times (in order to get more
% reliable estimates of the error).
for r=1:nr
    
    % Indices of the "unused" test samples.
    ind = [1:m];
    
    % Run through all the folds.
    for i=1:cvk
        
        % Indices of the current test data.
        tempind = randintex(ns(i),1,length(ind));
        testind = ind(tempind);
        trainind = indAll;
        trainind(testind) = [];
        ind(tempind) = [];
        
        % Infer the Boolean functions.
        [Fhat,OptE,Et] = bnBestFit(X(:,trainind),Y(:,trainind),w(trainind),k,F,X(:,testind),Y(:,testind),w(testind));
                
        Ehat = Ehat + Et;
        
    end % for i=1:cvk
    
    % Display something...
    disp([num2str(r),'/',num2str(nr)]);
    
end % for r=1:nr

% Normalize the error.
Ehat = Ehat/(nr*m);
