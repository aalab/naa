function Y = myunique(X,onlyUniqueElements)

% MYUNIQUE returns the unique values and their counts
% Copyright Sohan Seth sohan.seth@hiit.fi

if nargin == 1
    onlyUniqueElements = false;
end

[A, B, C] = unique(sort(X));
Y = sortrows([A(:), [B(2:end); length(C)+1] - B(:)],2);

if onlyUniqueElements
    Y = Y(Y(:,2) == 1);
end