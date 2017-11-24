time = linspace(-10,50,1000);
Ts = diff(time(1:2));
Fs = 1/Ts;
signal_1 = sinc(time);
signal_2 = sinc(time - 30)+0.05*randn(1,1000);

figure,
plot(time,signal_1);
figure,
plot(time,signal_2);

[acor,lag] = xcorr(signal_1,signal_2);
[~,I] = max(abs(acor));
lagDiff = lag(I);
timeDiff = lagDiff/Fs;

figure,
plot(time + timeDiff,signal_2)