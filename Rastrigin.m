function f1 = Rastrigin(x)

D  = numel(x);
f1 = sum(x.^2)+10*(D-sum(cos(2*pi*x)).^2);  