function [P,UX,UY] = compressible_flow_solver2D();
% 2D eliptic FD and FV solver
% compressible flow
% |---*---|---*---|---*---|
% Note : number of interface = number of cells + 1

% close all
% clear all

%% User input 

% disp('Welcome in 1D Fully Compressible Flow Simulator by RB Arbarim - 4573900')
% disp(' ')
% LX = str2double(input('Insert Reservoir Length (m) = ','s')); %length of the reservoir (m)
% LY = str2double(input('Insert Reservoir Width (m) = ','s')); %width of the reservoir (m)
% NX = round(str2double(input('Insert Number of Grid Cells in x direction = ','s')));
% NY = round(str2double(input('Insert Number of Grid Cells in y direction = ','s')));
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

LX = 10;
LY = 10; % LY becomes 1 in 1D
NX = 5;
NY = 5; % NY becomes 1 in 1D
t = 10e1;
nt = 10;
rho0 = 1000;
phi0 = 0.1;
P0 = 1;
cf = 1;
cr = 1;
well = 1;
nwell = 2;
grid = [1 NX*NY];
pi = [1000 1000];
pw = [1 0];
k = 10^-13; % absolute permeability in m2
visw = 10^-3; % water viscosity
homogenous = 1;

%% Initialization
dx = LX/NX; %grid size in x direction
x = linspace(dx/2,LX-dx/2,NX); %location of grid center
xi = linspace(0,LX,NX+1); %location of interfaces

dy = LY/NY; %grid size
y = linspace(dy/2,LY-dy/2,NY); %location of grid center
yi = linspace(0,LY,NY+1); %location of interfaces

dt = t/nt;
lambda = zeros(NY,NX);
lambdaharx = zeros(NY,NX+1);
lambdahary = zeros(NY+1,NX);
TX = zeros(NY,NX+1);
TY = zeros(NY+1,NX);
[xx yy] = meshgrid(x,y);

switch homogenous
    case(1)
    lambda = ones(NY,NX); % input lambda
    case(0)
    lambda = rand(NY,NX);
end

%% Time Loop

% Define initial condition
% P = ones(NX*NY,1)*P0;
P = zeros(NX*NY,1);
I = zeros(NY,NX);
[rho, drho] = computedensity(P,cf,rho0,P0);
[phi, dphi] = computeporosity(P,cr,phi0,P0);
rhom = vec2mat(rho,NX);;
[TX,TY] = computetransmisibility2D(rhom.*lambda,dx,NX,dy,NY);
[dummy_x,dummy_y,lambdaharx,lambdahary] = computetransmisibility2D(lambda,dx,NX,dy,NY);

q = computewellflux2D(pw,lambda,pi,grid,P);
uxcell = cell(nt,1);
uycell = cell(nt,1);

for Ndt = 1:round(nt)
% Newton Loop
converged = 0;
rhon = rho;
phin = phi;
    % to find P at n+1
    Rfs = computeresidual2D(rhon,phin,rho,phi,TX,TY,dt,P,q,NX,NY);
    while converged == 0
        % construct Jacobian matrix
        % J = A + C 
        A = zeros(NX*NY,NX*NY);
%         q = zeros(NX*NY,1);

      for i = 1 : NX
        for j = 1:NY
            gridindex(j,i) = (j-1)*NX+i;
            [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
        if i > 1 % there is left neighbor
            A(I,I)= A(I,I)+TX(j,i);
            A(I,IW)= -TX(j,i);
        end
        if i < NX % there is right neighbor
            A(I,I)= A(I,I)+TX(j,i+1);
            A(I,IE)= -TX(j,i+1);
        end
        if j > 1 % there is top neighbor
            A(I,I)= A(I,I)+TY(j,i);
            A(I,IN)= -TY(j,i);
        end
        if j < NY % there is bottom neighbor
            A(I,I)= A(I,I)+TY(j+1,i);
            A(I,IS)= -TY(j+1,i);
        end
        end
      end

    vec = (dphi.*rho + drho.*phi)/dt;
    C = diag(vec);
    
    % wells
    % Q = rho*lambda*PI*(pwell-p(i))
    % W(i,i) = lambda(i)*PI(i)*rho(i);
    W = zeros(NX*NY,NX*NY); % only non zero at well location
    
   for i = 1:length(pw);
       W(grid(i),grid(i))= W(grid(i),grid(i))+lambda(grid(i))*pi(i)*rho(grid(i));
   end
    
    Jfs = A+C+W;
    dp = Jfs\Rfs;
    
    %find pressure at nu+1
    P = P + dp;
    
    % update fluid and rock properties
    [rho, drho] = computedensity(P,cf,rho0,P0);
    [phi, dphi] = computeporosity(P,cr,phi0,P0);
    rhom = vec2mat(rho,NX);
    [TX,TY] = computetransmisibility2D(rhom.*lambda,dx,NX,dy,NY);
    q = computewellflux2D(pw,lambda,pi,grid,P);
    
    % compute residual
    Rfs = computeresidual2D(rhon,phin,rho,phi,TX,TY,dt,P,q,NX,NY);
     
        if norm(Rfs,2) < 1e-6;
            converged = 1;
        end
    end
    
[UX,UY] = computevelocity2D(lambdaharx,lambdahary,dx,dy,P,NX,NY,i,j);
p(Ndt,:) = P;
uxcell{Ndt} = UX;
uycell{Ndt} = UY;

figure(2)
hold on
surf(xx,yy,vec2mat(P,NX))
end
xlabel('x location')
ylabel('y location')
zlabel('Pressure')
title('Fully Compressible Pressure Solution')
hold off

%% Plot solution

P = vec2mat(P,NX);
% [xx yy] = meshgrid(x,y);
figure
surf(xx,yy,P)
xlabel('x grid location')
ylabel('y grid location')
zlabel('Pressure')
title('Compressible Pressure Solution')

figure
subplot(1,2,1)
surf(UX)
xlabel('x grid location')
ylabel('y grid location')
zlabel('Velocity')
xlim([0 NX])
ylim([0 NY])
title('Compressible Velocity Solution in x Direction')

subplot(1,2,2)
surf(UY)
xlabel('x location')
ylabel('y location')
zlabel('Velocity')
xlim([0 NX])
ylim([0 NY])
title('Compressible Velocity Solution in y Direction')

end