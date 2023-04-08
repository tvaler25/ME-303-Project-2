clear 
clc
close all

L= 1;         % x in (0,L)
T= 1;       % t in (0,T)
k=2;    % conductivity will not change
J = 10000;%iterations for analytical 

%Find Cn from n=1 to n=Jmax
Cn = zeros(J, 1);
Cn(1) = -4; %using limits since term evaluates to 0/0
for n=2:1:J
    Cn(n) = 2*((n^2)-((-1)^n)*(2-(3*(n^2))))/(n*((n^2)-1));
end
Cn = Cn/pi;

dt_array=      [1/1000, 1/2000, 1/5000, 1/10000, 1/20000, 1/50000];
dx_array=  [1/25, 1/20, 1/16, 1/12, 1/10, 1/5];

error = zeros(length(dt_array), length(dx_array)); %generate an array to stare error

for dx_number=1:length(dx_array)
    dx = dx_array(dx_number);
    N=round(L/dx);

    for dt_number=1:length(dt_array)

        dt = dt_array(dt_number)
        M=round(T/dt);

        F=k*dt/dx^2;
        if (1-2*F)>0 %check if unstable 
            
            %set matrix
            numer = zeros(N+1, M+1);
            exact = zeros(N+1, M+1);
    
            % Position of nodes
            x = linspace(0, L, N+1);
    
            % Initial Condition
            numer(:, 1) = cos(pi * x);
            numer(1,1)=0;
            numer(N+1,1)=2;
    
            % Explicit Scheme for Partial Difference Equation
            for j=1:M % time coordinate = j/M
                
                for i=2:N % space coordinate = i/N
                    numer(i, j+1) = numer(i, j) + F * (numer(i+1, j) - 2*numer(i, j) + numer(i-1, j));
                end
                
                numer(1, j+1) = 0; % BC left
                numer(N+1, j+1) = 2; % BC right: a time-varying one
            end
            
            for j=1:M+1 %run analytical solution 
    
                for i=1:N+1
                    %find exact temperature
                    exact(i, j) = exactTemp(((i-1)*dx), ((j-1)*dt), Cn);
                end
            end
            
            error(dt_number,dx_number) = rmse(exact, numer); %make error matrix
               
%             figure('Name','Numerical')
%             [X, Time] = meshgrid(0:dt:T, x);
%             surf(X,Time, numerical)
%             shading interp
%             colormap('jet')
%             xlabel('t'); 
%             ylabel('x'); 
%             zlabel('T(x,t)'); colorbar
        
        else
            disp("Unstable")
        end
    end
end

%% plot

figure('Name', 'Errors')
loglog(dt_array, error(:,1), 'LineWidth', 2, 'Marker', 'square');
hold on
for i=2:1:length(dx_array)
    plot(dt_array, error(:,i), 'LineWidth', 2, 'Marker', 'square');
end
grid;
hold off
title('Error In Numerical Solutions')
xlabel('Grid Spacing (dt)')
ylabel('RMS error')
legendText =  "dx= " + string(dx_array);
legend(legendText)
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