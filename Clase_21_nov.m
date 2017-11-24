
w = [3 7 5];
    
I = [3 2 4 1 3 8 4 0 3 8 0 7 7 7 1 2];


Mw = diag(w(1)*ones(1,numel(I)-1),-1)+diag(w(2)*ones(1,numel(I)),0)...
    + diag(w(3)*ones(1,numel(I)-1),+1);
p = double(Mw ~= 0);

corr_1 = conv(I,fliplr(w),'same');

corr_2 = sum(p.*(Mw - ones(16,1)*I).^2,2)';

s_II = sqrt(sum((p.*(ones(16,1)*I)).^2,2))';
s_ww = sqrt(sum(w.^2));
corr_3 = corr_1/s_ww./s_II

