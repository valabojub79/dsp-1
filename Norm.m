% X : Vector
% N : 1, 2, 3, ... , Inf

function Y = Norm(X,N)

if nargin < 2,
    N   = 2;    
    if nargin < 1,
        help Norm;
        return;
    end
end

switch N,
    case 1,
        Y = sum(abs(X));
    case 2,
        Y = sqrt(sum(X.^2));
    case inf,
        Y = max(abs(X));
    otherwise,
        Y = (sum(X.^N))^(1/N);
end
