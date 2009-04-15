function [ab, d] = bnAttractor(Avec)

% [ab,d] = bnAttractor(Avec) - attractors and distances to the attractors
% This function takes as input a state transition matrix (in vector form;
% see bnAsparse) and outputs two vectors, ab and d, which contain the
% attractors and distances to the attractors, respectively. In the vector
% ab, the attractors are numbered by negative numbers (-1, -2, -3, ...) and
% the corresponding basins are numbered by the corresponding positive
% numbers (1, 2, 3, ...). In the vector d, the distance from each state to
% its corresponding attractor is given. Naturally, if a state is on the
% attractor, its distance is 0. Both ab and d are of length 2^n, where n is
% the number of genes.
% Ilya & Harri; 09/11/02

nstates = length(Avec);
ab = zeros(1,nstates);
d = ab;
k = 1;                                          % k is the number of the attractor
ptr = 1;
abk = ab;
% wb = waitbar(0,'computing attractors and distances...');

while ptr <= nstates
%     waitbar(ptr/nstates,wb);
    %ptr = min(find(~ab));
    ab(ptr) = k;
    abkidx = 1;
    abk(abkidx) = ptr;
    dist = 1;
    d(ptr) = dist;                              % start with distance 1 and keep incrementing it
    nextstate = Avec(ptr);                      % makes one network step from ptr state
    
    while ~ab(nextstate)                        % as long as we haven't yet come back to a state marked by k
        ab(nextstate) = k;
        abkidx = abkidx + 1;
        abk(abkidx) = nextstate;
        dist = dist + 1;
        d(nextstate) = dist;
        nextstate = Avec(nextstate);
    end
    maxdist = d(nextstate) - 1;                 % maximum distance to that attractor
    
    if ab(nextstate) == k                       % we found a new attractor, so go through it again and change all states in it to -k
        while ab(nextstate) == k                % keep going until we run into an already marked state
            d(nextstate) = 0;                   % erase the distances because we're on the attractor
            ab(nextstate) = -k;                 % we mark the attractor with -k
            abkidx = abkidx - 1;
            nextstate = Avec(nextstate);        % keep going
        end
        %indbasin = find(ab == k);               % indices of current basin
        indbasin = abk(1:abkidx);
        %d(indbasin) = max(d(indbasin)) - d(indbasin) + 1;   % reverse the distances to the attractor
        d(indbasin) = maxdist - d(indbasin) + 1;   % reverse the distances to the attractor
        k = k + 1;
    elseif ab(nextstate) < 0                    % we found an old attractor
        %indbasin = find(ab == k);               % indices of current basin  (at this point, we still pretend we are going to a new attractor)
        indbasin = abk(1:abkidx);
        %d(indbasin) = max(d(indbasin)) - d(indbasin) + 1;   % reverse the distances to the attractor
        d(indbasin) = dist - d(indbasin) + 1;   % reverse the distances to the attractor
        ab(indbasin) = -ab(nextstate);     % we go back and change everything upstream to match the old attractor
    else                                        % we found an old basin
        %indbasin = find(ab == k);               % indices of current basin
        indbasin = abk(1:abkidx);
        %d(indbasin) = max(d(indbasin)) - d(indbasin) + 1 + d(nextstate);   % reverse the distances to the attractor and add the distance of the state in the old basin to which we arrived
        d(indbasin) = dist + d(nextstate) - d(indbasin) + 1;
        ab(indbasin) = ab(nextstate);           % we go back and change everything upstream to match the old basin
    end 
    while (ptr <= nstates) & ab(ptr)
        ptr = ptr + 1;
    end
end
% close(wb)