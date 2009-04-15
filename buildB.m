function B = buildB(n,mu)

% SYNTAX: B = buildB(n,mu)
% This function builds the transform matrix B(n,mu), which is used to test
% whether a Booealn function is A^mu.
% Ilya Shmulevich; Aug. 29, 2002

all_combs = nchoosek(2:2^n,mu)-1;                   % we don't want the zero-vectors here, so we start from 2
ncombs = length(all_combs);
conjs = zeros(ncombs,1);

for i = 1:ncombs,
    cjs = all_combs(i,1);
    for j = 2:mu,
        cjs = bitand(cjs,all_combs(i,j));
    end
    conjs(i) = cjs;
end

b = all_combs(logical(~conjs),:)+1;                 % add 1 for indexing
B = zeros(size(b,1),2^n);

for i = 1:size(B,1),
    B(i,b(i,:))=1;
end