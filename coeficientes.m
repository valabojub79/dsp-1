clc;
close all;
clear all;

N = 500; 

f_max = 180; 

%Ts = 1/1800; 
Ts = 1/2/f_max;

f_i = linspace(0,f_max,N);
 
w_i = 2*pi*Ts*f_i;

M_d_i = (f_i >= 45)&(f_i <= 55);

W_i = ones(1,N);

p = 2;

boundaries = ones(7,1)*[-10 10];

ObjFun = @(x) ObjFunIIR(x,w_i,M_d_i,W_i,p);

%[X_opt,F_opt] = PSO_IIR(ObjFun,boundaries)

[X_opt,F_opt] = UPSO(ObjFun,boundaries)



[~,M_i]=ObjFunIIR(X_opt,w_i,M_d_i,W_i,p);
%[~,M_i]=ObjFunIIR(X_opt); %nargin

plot(f_i,M_i,'r','LineWidth',1.5);
hold on;
plot(f_i,M_d_i,'k','LineWidth',1.5);
