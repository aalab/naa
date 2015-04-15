function val = myrands(clustA, clustB)

% val = myrands(clustA, clustB) computes the Rand's index between two 
% clustering assignments given in clustA and clustB
% Author: Sohan Seth, sohan.seth@hiit.fi

A = double(bsxfun(@eq, clustA, clustA')) -2 * (eye(length(clustA)));
B = double(bsxfun(@eq, clustB, clustB')) -2 * (eye(length(clustB)));
A = A(:); A(A == -1) = [];
B = B(:); B(B == -1) = [];

a = sum(A(:) .* B(:));
b = sum((1 - A(:)) .* (1 - B(:)));
%c = sum(A(:) .* (1 - B(:)));
%d = sum((1 - A(:)) .* B(:));

val = (a + b) / nchoosek(length(clustA), 2) / 2; %(a + b) / (a + b + c + d);