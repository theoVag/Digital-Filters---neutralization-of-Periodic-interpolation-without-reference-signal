%Digital filter assignments
%author: Theodoros- Panagiotis Vagenas
clear all
close all
n=200000;
M=100;
sigmav2=0.54;
%%generate white noise
v=rand(n,1);
v=sqrt(sigmav2)*v;
v=v-mean(v);
%%generate input signal
fo=1/4;
ph=pi/2;
A=4.2;
delta=10;
t=1:n;
x(t,1)=A*(sin(2*pi*fo*t+ph)+cos(4*pi*fo*t+ph)+cos(7*pi*t+ph/3));
s=x+v;
u = zeros(size(s));
u(delta + 1:end) = s(1:end - delta);
%% calculate correlations
[r,lags]=xcorr(u,u,M+1,'unbiased');
r=r(lags>=0);
R = toeplitz(r(1:end - 1));  % autocorrelation with r(0)
rbar = r(2:end);  % rbar without r(0).
rbar = rbar(:); % rbar column vector

%Wiener coefficients
%tic
wo=R\rbar; 
%toc
%%levinson function + compare with matlab
[a, G, L, Dp] = LevinsonDurbin_iterative(M, r);
[mat_a,mat_p,mat_g]=levinson(r,M);

fprintf('Filter coefficient norm(a-mat_a): %e \n', norm(a-mat_a.'));
fprintf('G norm(G-mat_g):  %e \n', norm(G-mat_g));
fprintf('Prediction Error power norm(Dp(end) - mat_p): %e \n', norm(Dp(end) - mat_p));

%%Gram-Schmidt orthogonalization algorithm
b=zeros(n,M+1);
b(1:M)=u(1:M);

for i=M+1:n
   b(i,:) = L * u(i:-1:i-M);
end
D=diag(Dp);
%%crosscorrelation of b,s
p=zeros(M+1,1);
for m=1:M+1
   [rb,lags]=xcorr(b(:,m),s,M,'unbiased');
   p(m) = rb(lags == 0);
end
%optimal param for joint process estimator
%tic
go=D\p;
%toc
%verify  L.''*go = wo (norm)
fprintf('norm(L.''*go-wo):%e\n',norm(L.'*go-wo));

%verify L.''*go = wo (mse)
fprintf('MSE L.''*go = wo :%e\n', immse(L.' * go, wo));

%wiener error
y=zeros(n,1);
for i=M+2:n
   y(i)=wo'*u(i-1:-1:i-M-1);
end
y = y(1:length(x));
e=s-y;
MSE1=(e-v).^2;
%Joint process estimator error
y=zeros(n,1);
y(1:M)=u(1:M);
for i=M+1:n
   y(i)=b(i,:)*go;
end
e=s-y;
MSE2=(e-v).^2;

figure(1)
subplot(1, 2, 1);
plot(MSE1)
xlabel('Time steps')
ylabel('Error');
title('Wiener MSE')

subplot(1, 2, 2);
plot(MSE2)
xlabel('Time steps')
ylabel('Error');
title('Joint Process Estimator MSE')
