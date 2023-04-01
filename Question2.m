
clear 
clc
close all

L= 1;         % x in (0,L)
T= 0.1;       % t in (0,T)
k=2;    % conductivity
N=20;   % cut space into N sections
M=2500; % cut time  into M sections
dx=L/N; dt=T/M; % grid spacing

F=k*dt/dx^2;


temp = zeros(N+1, M+1);

% Position of nodes
x = linspace(0, L, N+1);
 
% Initial Condition
temp(:, 1) = cos(pi * x);
 
% Explicit Scheme for Partial Difference Equation
for j=1:M % time coordinate = j/M
    
    for i=2:N % space coordinate = i/N
        temp(i, j+1) = temp(i, j) + F * (temp(i+1, j) - 2*temp(i, j) + temp(i-1, j));
    end
    temp(1, j+1) = 0; % DBC left
%     temp(1, j+1) = 2 + temp(2, j+1); % NBC left
    
%     temp(N+1, j+1) = 0; % DBC right: a constant BC
    temp(N+1, j+1) = 2; % DBC right: a time-varying one
end
 



%% plot
figure()
[T, X] = meshgrid(0:dt:T, x); 
surf(T, X, temp); 
shading interp
colormap('hot')
xlabel('t'); 
ylabel('x'); 
zlabel('T(x,t)'); colorbar


