% X : Vector
% N = 1 Norma 1
% N = 2 Norma 2
% N = 3 Norma MAX

function Y = Norm(X,N)
if (N == 1)
    Y = sum(abs(X));
elseif (N == 2)
    Y = sqrt(sum(X.^2));
elseif (N == 3)
    Y = max(abs(X));
else 
    disp(' Error ')
end
