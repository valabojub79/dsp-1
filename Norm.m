% X : Vector
% N : 1, 2, 3, ... , Inf

function Y = Norm(X,N)

switch N,
    case 1,
        Y = sum(abs(X));
    case 2,
        Y = sqrt(sum(X.^2));
    case inf,
        Y = max(abs(X));
    case -inf
        Y = min(abs(X));
    otherwise,
        Y = (sum(X.^N))^(1/N);
end