
clear 
clc
close all

L= 1;         % x in (0,L)
T= 0.1;       % t in (0,T)
k=2;    % conductivity
N=20;   % cut space into N sections
M=2500; % cut time  into M sections
J=100; % use J iterations of the summation
dx=L/N; dt=T/M; % grid spacing

F=k*dt/dx^2;

%Find Cn from n=1 to n=J
Cn = zeros(J, 1);
Cn(1) = -4; %using limits since term evaluates to 0/0
for n=2:1:J
    Cn(n) = 2*((n^2)-((-1)^n)*(2-(3*(n^2))))/(n*((n^2)-1));
end
Cn = Cn/pi;

temp = zeros(N+1, M+1);
exact = zeros(N+1, M+1);

% Position of nodes
x = linspace(0, L, N+1);
 
% Initial Condition
temp(:, 1) = cos(pi * x);
 
% Explicit Scheme for Partial Difference Equation
for j=1:M % time coordinate = j/M
    
    for i=2:N % space coordinate = i/N
        temp(i, j+1) = temp(i, j) + F * (temp(i+1, j) - 2*temp(i, j) + temp(i-1, j));
        
        %find exact temperature
        exact(i, j) = exactTemp((i*dx), (j*dt), Cn);
    end
    temp(1, j+1) = 0; % DBC left
%     temp(1, j+1) = 2 + temp(2, j+1); % NBC left
    
%     temp(N+1, j+1) = 0; % DBC right: a constant BC
    temp(N+1, j+1) = 2; % DBC right: a time-varying one
end


%% plot
figure(1)
[T, X] = meshgrid(0:dt:T, x); 
surf(T, X, temp); 
shading interp
colormap('hot')
xlabel('t'); 
ylabel('x'); 
zlabel('T(x,t)'); colorbar


function exactTemp = exactTemp (x, t, Cn)
    sum = 0;
    for n=1:length(Cn)
        termN = Cn(n)*sin(n*pi*x)*exp(-2*(n^2)*(pi^2)*t);
        sum = sum + termN;
    end
    exactTemp = sum + (2*x);
end