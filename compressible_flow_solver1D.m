% function [P,U] = compressible_flow_solver1D();
% 1D eliptic FD and FV solver
% compressible flow
% |---*---|---*---|---*---|
% Note : number of interface = number of cells + 1

close all
clear all
clc
format short

%% Input Data
% disp('Welcome in 1D Fully Compressible Flow Simulator by RB Arbarim - 4573900')
% disp(' ')
% L = str2double(input('Insert Reservoir Length (m) = ','s')); %length of the reservoir (m)
% N = round(str2double(input('Insert Number of Grid Cells = ','s')));
% t = str2double(input('Insert time = ','s')); % length of time 
% nt = str2double(input('Insert time grid = ','s')); % time grid
% rho0 = str2double(input('Insert initial density = ','s')); 
% phi0 = str2double(input('Insert initial porosity = ','s')); 
% P0 = str2double(input('Insert initial pressure = ','s')); 
% cf = str2double(input('Insert fluid compressibility = ','s')); 
% cr = str2double(input('Insert rock compresibility = ','s')); 
% % well = str2double(input('Do you want to put wells (Press 1 for Yes, press 0 for No) = ','s'));
% % 
% % if well == 1
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
%    
%    homogenous = str2double(input('Press 1 for homogenous reservoir, 0 for heterogenous reservoir = ','s')); % case homogenous or heterogenous

% solution = str2double(input('Press 1 for implicite solution, 0 for explicite solution = ','s')); % case homogenous or heterogenous

%% Default input
L = 100;
N = 100;
t = 1e5;
nt = 10;
rho0 = 1000;
phi0 = 0.1;
P0 = 2e8;
cf = 2e-8;
cr = 1e-8;
k = 10^-13; % absolute permeability in m2
visw = 10^-3; % water viscosity
densw = 1000; % water density in kg/m3
nwell = 2;
grid = [1 N];
pi = [1e5 1e5];
pw = [3e8 1e8];
homogenous = 1;

%% Initialization
dx = L/N; %grid size
dt = t/nt;
x = linspace(dx/2,L-dx/2,N); %location of grid center
xi = linspace(0,L,N+1); %location of interfaces

lambda = zeros(N,1);
lambdahar = zeros(N+1,1);
T = zeros(N+1,1);

switch homogenous
    case(1)
    lambda(1:N)= k/visw; % input lambda
    case(0)
    lambda(1:N) = k/visw*rand(N,1);
end

%% Time Loop

% Define initial condition
P = P0*ones(N,1);
[rho, drho] = computedensity(P,cf,rho0,P0);
[phi, dphi] = computeporosity(P,cr,phi0,P0);
T = computetransmisibility1D(rho.*lambda,dx,N);
[dummy,lambdahar] = computetransmisibility1D(lambda,dx,N);
time = cell(nt,1);
[q,dqdP] = computewellflux1D(pw,lambda,pi,grid,P);
display('start');

for Ndt = 1:round(nt)
% Newton Loop
converged = 0;
rhon = rho;
phin = phi;
    % to find P at n+1
%     R = computeresidual1D(rhon,phin,rho,phi,T,dt,P,q,N);
    [R,J] = computeresidual1D(rhon,phin,rho,drho,phi,dphi,T,dt,P,q,lambda,pi,N,grid);

    while converged == 0
    dP = J\-R;
    
    %find pressure at nu+1
    P = P + dP;
    
    % update fluid and rock properties
    [rho, drho] = computedensity(P,cf,rho0,P0);
    [phi, dphi] = computeporosity(P,cr,phi0,P0);
    T = computetransmisibility1D(rho.*lambda,dx,N);
    [q,dqdP] = computewellflux1D(pw,lambda,pi,grid,P);
    
    % compute residual
    [R,J] = computeresidual1D(rhon,phin,rho,drho,phi,dphi,T,dt,P,q,lambda,pi,N,grid);
    
        if norm(R,2) < 1e-6;
            converged = 1;
        end
    end
    
U = computevelocity1D(lambdahar,dx,P,N);
p(Ndt,:) = P;
u(Ndt,:) = U;
time{Ndt}= strcat('time step : ',num2str(Ndt));

hold on
figure(1)
plot(x,P)
end
title('Fully Compressible Flow Solver')
xlabel('x Location')
ylabel('Pressure')
% legend(time)
hold off

% for Ndt = 1 : round(nt)
%     figure(2)
%     hold on
%     plot(xi,u(Ndt,:))
% end
% xlabel('x location')
% ylabel('Velocity')
% title('Fully Compressible Flow Solver')
% legend(time)
% hold off
% 
% %% additional plot
% Ndt = [1:nt];
% timeindex = [Ndt(1) round(0.25*Ndt(end)) round(0.5*Ndt(end)) round(0.75*Ndt(end)) Ndt(end)];
% timelegend = cell(length(timeindex),1);
% 
% for i = 1 : length(timeindex)
%     timelegend{i}= strcat('time step ',num2str(timeindex(i)));
%     figure(3)
%     hold on
%     plot(x,p(timeindex(i),:))
% end
% legend(timelegend)
% hold off
% xlabel('x location')
% ylabel('Pressure')
% title('Fully Compressible Flow Solver')
% 
% for i = 1 : length(timeindex)
%     timelegend{i}= strcat('time step ',num2str(timeindex(i)));
%     figure(4)
%     hold on
%     plot(xi,u(timeindex(i),:))
% end
% legend(timelegend)
% hold off
% xlabel('x location')
% ylabel('Velocity')
% title('Fully Compressible Flow Solver')
% 
% % return