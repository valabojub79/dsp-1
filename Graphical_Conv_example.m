% Example from p. 36 of Lecture Notes
t = [-4 -3 -2 -1 0 1 2 3 4];
x = [ 0 0 1 1 1 1 1 0 0];
h = [ 0 0 0 0 1 1 1 0 0];
y = [ 0 0 0 0 0 0 0 0 0];

yc = 1;

for n=min(t):max(t)
    pause(3);
    
    % flip h
    ht = fliplr(h);
    
    if n<0
        % shift to the left
        ht = [ht(-n+1:length(h)) zeros(1,-n)];
    else
        % shift to the right
        ht = [zeros(1,n) ht(1:length(h)-n)];
    end
    
    y(yc) = sum(x.*ht);
    yc = yc + 1;
    
    subplot(2,1,1);
    stem(t,x,'linewidth',2);
    hold on;
    stem(t,ht,'filled','k','linewidth',2);
    hold off;
    xlabel('t');
    legend('x[n]','h[n-k]');
    title(['n= ' num2str(n)]);
    
    subplot(2,1,2);
    stem(t,y,'linewidth',2);
    xlabel('t'); ylabel('y[n]');
end