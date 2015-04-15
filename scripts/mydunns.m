function val = mydunns(d, clustInd)

% val = mydunns(d, clustInd) computes dunn's index given distance matrix d
% and cluster assignments clustInd
% Author: Sohan Seth, sohan.seth@hiit.fi

% for easy debugging
% [clustInd, ind] = sort(clustInd);
% d = d(ind, ind);

I = double(bsxfun(@eq, clustInd, clustInd'));
d_intra = d .* I; 
d_inter = d .* (1 - I) + max(d(:)) * I; % later part added to avoid zeros as minimum
val = min(d_inter(:))/max(d_intra(:));