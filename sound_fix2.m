%Digital filter assignments
%author: Theodoros- Panagiotis Vagenas
clear all
close all

M=100;
delta=100;
load music.mat
n=length(s);
d=s;
u=zeros(size(d));
u(delta+1:end)=s(1:end-delta);
[r,lags]=xcorr(u,u,M+1,'unbiased');
r=r(lags>=0);
R=toeplitz(r(1:end-1));
rbar=r(2:end);
wo=R\rbar;
PM=r(1)-rbar.'*wo;

[a, G, L, Dp] = LevinsonDurbin_iterative(M, r);
p=zeros(M+1,1);
for i=1:M+1
    b=filter(L(i,1:i),1,u);%function filter to reduce memory allocation
    b(1:i)=u(1:i);
    [rb,lags]=xcorr(b,s,1,'unbiased');
    p(i)=rb(lags==0);
end
D=diag(Dp);
go=D\p;

y=zeros(n,1);
for i=1:M+1
    b=filter(L(i,1:i),1,u);
    b(1:i)=u(1:i);
    y=y+go(i)*b;
end
e = s - y;
sound(e, fs);
