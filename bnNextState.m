function y = bnNextState_new(x,F,varF,nv)

% y = bnNextState(x,F,varF,nv) - one step of a Boolean Network
% This function performs one step of a BN. The current state is x (binary
% vector of 1 x n, where n is # of genes). The other parameters are the
% same as in bnA. Note: perturbation is not supported here.

% Ilya Shmulevich; Aug. 28, 2002
% Modified May 12, 2003 by HL.

n = length(x); % number of genes
y = zeros(size(x)); % initialize next state

b = 2.^[size(varF,1)-1:-1:0]';
for i=1:n
    if nv(i) ~= 0
        y(i) = F(x(varF(1:nv(i),i))*b(end-nv(i)+1:end)+1,i);
    else
        y(i) = x(i);
    end
end
