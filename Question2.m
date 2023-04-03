
clear 
clc
close all

L= 1;         % x in (0,L)
T= 10;       % t in (0,T)
k=2;    % conductivity

%grid spacing
%dt=0.03;

N=50;   % cut space into N sections
M=150000; % cut time  into M sections
J=1000; % use J iterations of the summation
dx=L/N; 
dt=T/M; % grid spacing

F=k*dt/dx^2;

if (1-2*F)>0

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
    temp(1,1)=0;
    temp(N+1,1)=2;
     
    % Explicit Scheme for Partial Difference Equation
    for j=1:M % time coordinate = j/M
        
        for i=2:N % space coordinate = i/N
            temp(i, j+1) = temp(i, j) + F * (temp(i+1, j) - 2*temp(i, j) + temp(i-1, j));
            
            
        end
        temp(1, j+1) = 0; % DBC left
        
        temp(N+1, j+1) = 2; % DBC right: a time-varying one
    end
    
    for j=1:M+1
    
        for i=1:N+1
            %find exact temperature
            exact(i, j) = exactTemp(((i-1)*dx), ((j-1)*dt), Cn);
        end
    end
    
    
    %% plot
    
    figure('Name','Analytical')
    [X, Time] = meshgrid(0:dt:T, x);
    surf(X,Time, exact)
    shading interp
    colormap('jet')
    xlabel('t'); 
    ylabel('x'); 
    zlabel('T(x,t)'); colorbar
    
    figure(2)
    [Time, X] = meshgrid(0:dt:T, x); 
    surf(Time, X, temp); 
    shading interp
    colormap('jet')
    xlabel('t'); 
    ylabel('x'); 
    zlabel('T(x,t)'); colorbar
    
    
   
    sliced_time_indices=round([0.001,0.01,0.1,10]/dt);
    num_slices=length(sliced_time_indices);
    
    slices=zeros(N+1,num_slices*2);
    
    for slice_number=1:num_slices %creating slices to compare profiles of ndiffernt time values
        slices(:,slice_number)=temp(:, sliced_time_indices(slice_number));
        slices(:,slice_number+4)=exact(:, sliced_time_indices(slice_number));
    end

    position_x=linspace(0,L,N+1)'


    
    %plots position
    figure('Name','2D Comparison')
    plot(position_x,slices(:, 1),'LineWidth', 2);
    hold on
    for i=2:1:num_slices*2
        plot(position_x, slices(:, i), 'LineWidth', 2);
    end
    hold off
    grid
%     legendText =  "F: " + string(frontStiffnessDisplay)+ ", R: " + string(rearStiffnessDisplay) + " (kN/rad)";
    %legend(legendText)
    xlabel('Position X')
    ylabel('Temp')

    
    
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

  