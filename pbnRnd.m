function [F,varF,cij] = pbnRnd(n,nf,nv)

% [F,varF,cij] = pbnRnd(n,nf,nv) - random BN/PBN
% This function creates a random PBN with n nodes. The number of functions
% per node is defined in nf. Each function may have different number of
% variables, defined in nv. Function returns all the necessary parameters
% (excluding nf and nv), to define an independent PBN.
% INPUT:
% n     - The number of nodes in a random pbn.
% nf    - The number of functions per node. nf is a 1-by-n vector where
%         each element of nf is an integer larger than 0. If nf is a vector
%         of ones, then this function generates a random Boolean network.
% nv    - The number of variables in each Boolean function. The length of
%         nv must be equal to the sum of the elements in vector nf. Let cnf
%         be the cumulative sum of the number of functions in the network
%         (starting from the first node), i.e., cnf = cumsum(nf);. Then,
%         the number of variables of the functions for the first node are
%         nv(1:cnf(1)), for the second node nv(cnf(1)+1:cnf(2)), ..., and
%         for the last node nv(cnf(n-1)+1:cnf(n)).
% OUTPUT:
% F     - The functions for the random network. F has size
%         2^max(nv)-by-sum(nf). Truth tables of the functions for the first
%         node are defined in F(:,1:cnf(1)), for the second node
%         F(:,cnf(1)+1:cnf(2)), ..., and for the last node
%         F(:,cnf(n-1)+1:cnf(n)). Since the length of the truth tables may
%         vary between different functions, only the first 2^nv(i) bits are
%         relevant in the i:th column of F. Remaining "unnecessary" bits in
%         each column are set to -1. Let the f = F(;,i) be the i:th
%         function in F, and assume that it is a function of three
%         variables xi, xj, and xk (variables are defined in varF(1:3,i).
%         Then, f(1) defines the output for the input vector (000).
%         Correspondinly, f(2) defines the output for the input vector
%         (001), where xi = xj = 0 and xk = 1, i.e., where the third input
%         variable is equal to one. As another example, f(6) defines the
%         output for the input vector (101), where xi = 1, xj = 0, and
%         xk = 1.
% varF  - The variables of each function. varF is a max(nv)-by-sum(nf)
%         matrix. Variables of the functions for the first node are
%         varF(:,1:cnf(1)), for the second node varF(:,cnf(1)+1:cnf(1)),
%         ..., and for the last node varF(:,cnf(n-1)+1:cnf(n)). Since the
%         number of variables may vary between different functions, only
%         the first nv(i) elements are relevant in the i:th column of varF.
%         Remaining "unnecessary" elements in each column are set to -1.
% cij   - The selection probabilities for the functions. cij is a
%         max(nf)-by-n matrix. The selection probabilities of the functions
%         for the i:th node are cij(:,1:i). Since the number of functions
%         may vary between nodes, only the first nf(i) elements are
%         relevant in the i:th column of cij. Remaining "unnecessary"
%         elements in each column are set to -1.
%
% EXAMPLES:
%
% Generate a random BN with 10 nodes, each having a 4-variable Boolean
% function.
% n = 10;
% nf = ones(1,n);
% nv = 4*ones(1,n);
% [F,varF,cij] = pbnRnd(n,nf,nv);
%
% Generate a random PBN with 10 nodes, each node having 2 functions, and
% the number of variables in each function being equal to 3.
% n = 10;
% nf = 2*ones(1,n);
% nv = 3*ones(1,sum(nf));
% [F,varF,cij] = pbnRnd(n,nf,nv);
%
% Generate a 10 node BN with ("truncated") scale-free topology (parameter
% gamma is equal to 3).
% n = 10;
% nf = ones(1,n);
% p = [1:n].^(-3);
% p = p/sum(p);
% nv = zeros(1,n);
% for i=1:n, nv(i) = pmfrnd(p);, end
% [F,varF,cij] = pbnRnd(n,nf,nv);
%
% Generate a random PBN with 10 nodes, each node having 2, 3, or 4
% functions (selected randomly from a pdf), and the number of variables in
% each function being 2, 3 or 4 (selected randomly from another pdf).
% n = 10;
% p1 = [0.3 0.4 0.3];
% p2 = [0.5 0.25 0.25];
% nf = zeros(1,n);
% for i=1:n, nf(i) = pmfrnd(p1) + 1; ,end
% nv = zeros(1,sum(nf));
% for i=1:sum(nf), nv(i) = pmfrnd(p1) + 1; ,end
% [F,varF,cij] = pbnRnd(n,nf,nv);

% Functions used: randintex.m

% Ilya Shmulevich; Aug. 14, 2001; Mar. 1, 2002;
% Modified: April 30, 2003 and May 14, 2003 by HL.


F = -ones(2^max(nv),sum(nf));
varF = -ones(max(nv),sum(nf));
cij = -ones(max(nf),n);

% Select the variables for each function and modify F.
for i=1:sum(nf)
    varF(1:nv(i),i) = sort(randintex(nv(i),1,n))';
    F(1:2^(nv(i)),i) = rand(2^nv(i),1)>0.5;
end % for i=1:sum(nf)

% Select the selection probabilities for each node.
for i=1:n
    cij(1:nf(i),i) = rand(nf(i),1);
    cij(1:nf(i),i) = cij(1:nf(i),i)/sum(cij(1:nf(i),i));
end % for i=1:n

