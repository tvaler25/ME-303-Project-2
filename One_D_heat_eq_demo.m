
clear 
clc
close all

L= 1;         % x in (0,L)
T= 0.5;       % t in (0,T)
k=10;    % conductivity
N=20;   % cut space into N sections
M=5000; % cut time  into M sections
dx=L/N; dt=T/M; % grid spacing
F=k*dt/dx^2;

temp = zeros(N+1, M+1);

% Position of nodes
x = linspace(0, L, N+1);
 
% Initial Condition
temp(:, 1) = sin(pi * x);
 
% Explicit Scheme for Partial Difference Equation
for j=1:M % time coordinate = j/M
    
    for i=2:N % space coordinate = i/N
        temp(i, j+1) = temp(i, j) + F * (temp(i+1, j) - 2*temp(i, j) + temp(i-1, j));
    end
    temp(1, j+1) = 0; % DBC left
%     temp(1, j+1) = 2 + temp(2, j+1); % NBC left
    
%     temp(N+1, j+1) = 0; % DBC right: a constant BC
    temp(N+1, j+1) = sin(pi * (j+1)/M); % DBC right: a time-varying one
end
 

%% animation
h = figure();
ax = gca;
ax.NextPlot = 'replaceChildren';

every_nth_frame = 50; %index freqeuncy for animation
num_frames = floor(M/every_nth_frame); %total number of frames

ani_temp = temp(:, 1:every_nth_frame:end); %temperature sampled at specified frequency
ax.XLim = [0, L];
ax.YLim = [min(ani_temp, [], "all"), max(ani_temp, [], 'all')];
xlabel('x'); ylabel('T(x,t)');

for idx = 1:num_frames
    plot(x, ani_temp(:, idx), '-r')
    drawnow
    frame = getframe(h);
    im{idx} = frame2im(frame);
end
close;

%save the gif, uncomment to save
% filename = "testAnimated.gif"; % Specify the output file name
% for idx = 1:num_frames
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0);
%     else
%         imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0);
%     end
% end


%% plot
figure()
[T, X] = meshgrid(0:dt:T, x); 
surf(T, X, temp); 
shading interp
colormap('hot')
xlabel('t'); ylabel('x'); zlabel('T(x,t)'); colorbar


