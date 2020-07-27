function [P,U] = slightly_compressible_flow_solver1D()
% 1D eliptic FD and FV solver
% incompresisble flow
% |---*---|---*---|---*---|
% Note : number of interface = number of cells + 1

% close all
% clear all
% format short

%% input data 
% 
% disp('Welcome in 1D Slightly Compressible Flow Simulator by RB Arbarim - 4573900')
% disp(' ')
% L = str2double(input('Insert Reservoir Length (m) = ','s')); %length of the reservoir (m)
% N = round(str2double(input('Insert Number of Grid Cells = ','s')));
% PL = str2double(input('Insert Pressure at Left Boundary (Pa) = ','s')); % Pressure at left boundary
% PR = str2double(input('Insert Pressure at Right Boundary (Pa) = ','s')); % Pressure at right boundary 
% t = str2double(input('Insert time = ','s')); % Pressure at right boundary
% nt = str2double(input('Insert time grid = ','s')); % Pressure at right boundary
% well = str2double(input('Do you want to put wells (Press 1 for Yes, press 0 for No) = ','s'));
% 
% if well == 1
%    nwell = str2double(input('How many wells do you want to put = ','s')); % Number of wells
%    fprintf('You have %g number of grid cell \n',N)
%    fprintf('In which grid cell you want to put wells (from 1 to %g ) = \n', N)
%    for i = 1 : nwell
%     fprintf('Well(%s) : \n',num2str(i))
%     grid(i) = round(str2double(input('Enter well location = ','s')));
%     if grid(i) > N
%         while grid(i) > N
%         disp('Invalid, well location exceeds number of grid cell')
%         grid(i) = str2double(input('Enter new well location = ','s'));
%         end
%     end
%     pi(i) = str2double(input('Enter productivity index = ','s'));
%     pw(i) = str2double(input('Enter wellbore pressure = ','s'));
%     disp('')
%    end
% end
% 
% homogenous = str2double(input('Press 1 for homogenous reservoir, 0 for heterogenous reservoir = ','s')); % case homogenous or heterogenous
% solution = str2double(input('Press 1 for implicite solution, 0 for explicite solution = ','s')); % case homogenous or heterogenous

%% Default input

N = 100; %number of grid cells
L = 1; %length of the reservoir (m)
PL = 1; %input boundary condition
PR = 0; %input boundary condition
t = 10e-2;
nt = 10;
well = 1;
nwell = 2;
pi = [1000 1000];
pw = [3 0];
grid = [1 N];
homogenous = 1;
solution = 0;

%% Initialization
% homogenous = 0;
dx = L/N; % space grid size
dt = t/nt; % time grid size 
x = linspace(dx/2,L-dx/2,N); %location of grid center
xi = linspace(0,L,N+1); %location of interfaces
ceff = 1; % effective compressibility in 1/Pa
k = 10^-13; % absolute permeability in m2
por = 0.1;
visw = 10^-3; % water viscosity
densw = 1000; % water density in kg/m3
D = k/(visw*por*ceff);
sigma = sqrt(2*D*t);

lambda = zeros(N,1);
lambdahar = zeros(N+1,1);
T = zeros(N+1,1);

switch homogenous
    case(1)
    lambda(1:N)= 1; % input lambda
    case(0)
    lambda(1:N) = rand(N,1);
end

[T,lambdahar] = computetransmisibility1D(lambda,dx,N);

%% Time Loop
% t = 1;
% dt = 0.001;
P = zeros(N,1);
time = cell(nt,1);

for Ndt = 1 : round(nt)
A = zeros(N,N);
q = zeros(N,1);

    for i = 1:N
        if i > 1 % there is left neighbor
            % T(i) * (P(i)-P(i-1))
            A(i,i)= T(i);
            A(i,i-1)= -T(i);
        end
        if i < N % there is right neighbor
            A(i,i)= A(i,i)+T(i+1);
            A(i,i+1)= -T(i+1);
        end
    end

switch well
    case(0)
    i = 1; % left boundary
    % T(1) *(P(1)-PL)
    A(i,i) = A(i,i)+ T(i);
    q(i) = q(i)+T(i)*PL;

    i = N; % right boundary
    % T(N+1) *(P(N)-PR)
    A(i,i) = A(i,i)+ T(i+1);
    q(i) = q(i)+T(i+1)*PR;
    
    case(1)
    [A,q] = wells1D(pi,pw,lambda,A,q,grid);
end

C = eye(N)*(ceff*por/dt);

% solution=0; %0 for explicit, 1 for implicite
switch solution
    case(1)
    P = C\(q-A*P+C*P);
    case(0)
    P = (C+A)\(q+C*P);
end

U = computevelocity1D(lambdahar,dx,P,N);
p(Ndt,:) = P;
u(Ndt,:) = U;
time{Ndt}= strcat('time step : ',num2str(Ndt));

figure(1)
hold on
plot(x,P)
end
xlabel('x location')
ylabel('Pressure')
legend(time)
title('Slightly Compressible Flow Solver')
hold off

for Ndt = 1 : nt
    figure(2)
    hold on
    plot(xi,u(Ndt,:))
end
xlabel('x location')
ylabel('Velocity')
title('Slightly Compressible Flow Solver')
legend(time)
hold off

%% additional plot
Ndt = [1:nt];
timeindex = [Ndt(1) round(0.25*Ndt(end)) round(0.5*Ndt(end)) round(0.75*Ndt(end)) Ndt(end)];
timelegend = cell(length(timeindex),1);

for i = 1 : length(timeindex)
    timelegend{i}= strcat('time step ',num2str(timeindex(i)));
    figure(3)
    hold on
    plot(x,p(timeindex(i),:))
end
legend(timelegend)
hold off
xlabel('x location')
ylabel('Pressure')
title('Slightly Compressible Flow Solver')

for i = 1 : length(timeindex)
    timelegend{i}= strcat('time step ',num2str(timeindex(i)));
    figure(4)
    hold on
    plot(xi,u(timeindex(i),:))
end
legend(timelegend)
hold off
xlabel('x location')
ylabel('Velocity')
title('Slightly Compressible Flow Solver')

end