function vfind = pbnFindFunc(F,varF,nf,nv,cij,p,states)

% vfind = pbnFindFunc(F,varF,nf,nv,cij,p,states) - alter the functions of a PBN
% This function takes as input a PBN, given by parameters F, varF, nf, nv,
% cij, and p. It alters every possible predictor (i.e. all predictors for
% all genes) by trying every possible Boolean function. So, it produces an
% "altered" PBN that is different from the given PBN by only one Boolean
% function. For each such altered PBN, it computes a steady-state
% distribution vector. The argument 'states' is a vector that specifies
% which of steady-state probabilities we're interested in seeing. For
% example, if states = [0 7], then we're interested in the stead-state
% prob. of (000) and (111). The output of this function is vfind, which is
% a matrix whose every row contains the steady-state probs. of the states
% given in 'states' and the last two columns contain the index of the
% altered predictor (column number in F) and the decimal value of the truth
% table of the (new) altered predictor.

% Function used: pbnA.m

% Ilya Shmulevich 01-31-02
% Modified May 14, 2003 by HL.

[nn,N]=size(F); % N is the number of all predictor functions in the given PBN.
n = log2(nn);   % number of variables
Fnew = F;
vfind = zeros(sum(2.^nv),2+length(states));
k = 1;
for i = 1:N,
    fprintf(1,'Modifying column %d/%d\n', i, N);
    for j = 1:2^(2^nv(i))
        f = bitget(j-1,[2^nv(i):-1:1]);
        Fnew(1:2^nv(i),i) = f';
        [A,v] = pbnA(Fnew, varF, nf, nv, cij, p);
        vfind(k,:) = [v(states+1)' i j-1];
        k = k + 1;
        Fnew = F;
    end
end
