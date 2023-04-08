format long
clear 
clc
close all

L= 1;         % x in (0,L)
T= 1;       % t in (0,T)
k=2;    % conductivity will not change




N=15;   % cut space into N sections
M=10000; % cut time  into M sections
iter=[1001, 1, 3, 6, 10, 18, 30, 60, 100, 180, 300, 600, 1000]; % use J iterations of the summation
dx=L/N; 
dt=T/M; % grid spacing
F=k*dt/dx^2;

if (1-2*F)>0
    
    error = zeros(length(iter) - 1, 1);
    
    for p=1:length(iter)
        
        J = iter(p);
        
        %Find Cn from n=1 to n=Jmax
        Cn = zeros(J, 1);
        Cn(1) = -4; %using limits since term evaluates to 0/0
        for n=2:1:J
            Cn(n) = 2*((n^2)-((-1)^n)*(2-(3*(n^2))))/(n*((n^2)-1));
        end
        Cn = Cn/pi;
        
        numerical = zeros(N+1, M+1);
        analytical = zeros(N+1, M+1);

        % Position of nodes
        x = linspace(0, L, N+1);

        % Initial Condition
        numerical(:, 1) = cos(pi * x);
        numerical(1,1)=0;
        numerical(N+1,1)=2;

        % Explicit Scheme for Partial Difference Equation
        for j=1:M % time coordinate = j/M

            for i=2:N % space coordinate = i/N
                numerical(i, j+1) = numerical(i, j) + F * (numerical(i+1, j) - 2*numerical(i, j) + numerical(i-1, j));


            end
            numerical(1, j+1) = 0; % DBC left

            numerical(N+1, j+1) = 2; % DBC right: a time-varying one
        end

        for j=1:M+1

            for i=1:N+1
                %find exact temperature
                analytical(i, j) = exactTemp(((i-1)*dx), ((j-1)*dt), Cn);
            end
        end
        if p==1
            exact = analytical;
        else
            error(p-1) = rmse(exact, analytical);
        end
        
        figure('Name','Analytical')
        [X, Time] = meshgrid(0:dt:T, x);
        surf(X,Time, numerical)
        shading interp
        colormap('jet')
        xlabel('t'); 
        ylabel('x'); 
        zlabel('T(x,t)'); colorbar
    end
    %% plot
iter(1) = [];
figure('Name', 'Errors')
loglog(iter, error, 'LineWidth', 2, 'color', 'red', 'Marker', 'square');
grid
title('error in numerical solutions')
xlabel('# of iterations in numerical solution')
ylabel('RMS error)')
    
else
    disp("Unstable")
end

function exactTemp = exactTemp (x, t, Cn)
    sum = 0;
    
    for n=1:length(Cn)
        termN = Cn(n)*sin(n*pi*x)*exp(-2*(n^2)*(pi^2)*t);
        sum = sum + termN;
    end
    exactTemp = sum + (2*x);
end

function rmse = rmse (A, B)
    diff = A - B;
    size = numel(diff);
    diffsqrd = diff.^2;
    sumdiffsqrd = sum(sum(diffsqrd));
    rmse = sqrt(sumdiffsqrd/size);
end