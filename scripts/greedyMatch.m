function indMatch = greedyMatch(d)

% indMatch = greedyMatch(d) matches row and column indices 
%   in a greedy fashion given dissimilarity matrix d. 
%   greedyMatch matches rows and columns with minimum distance
%   first and then remove this row/column from comparison
%
% Author: Sohan Seth, sohan.seth@hiit.fi

indMatch = zeros(1, size(d,2));
for count = 1:length(d)
   [values, indices] = min(d);
   [~, ind] = min(values);
   indMatch(ind) = indices(ind);
   d(indices(ind), :) = Inf;
   d(:, ind) = Inf;
end