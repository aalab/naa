function val = entDir(alpha)
% ENTDIR computes the entropy of Dirichlet distribution
%   VAL = ENTDIR(ALPHA) computes the entropy of each column ALPHA
%   and returns the values in vector VAL
%   
%   Copyright Sohan Seth sohan.seth@hiit.fi

val = sum(gammaln(alpha)) - gammaln(sum(alpha)) ...
    + (sum(alpha) - size(alpha, 1)) .* psi(sum(alpha)) ...
    - sum((alpha-1) .* psi(alpha));