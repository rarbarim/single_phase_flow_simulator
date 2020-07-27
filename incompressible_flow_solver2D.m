function [P,UX,UY] = incompressible_flow_solver2D()
% 2D eliptic FD and FV solver
% incompresisble flow
% |---*---|---*---|---*---|
% Note : number of interface = number of cells + 1

% close all
% clear all
% format short
%% input data 

% disp('Welcome in 1D Incompressible Flow Simulator by RB Arbarim - 4573900')
% disp(' ')
% LX = str2double(input('Insert Reservoir Length (m) = ','s')); %length of the reservoir (m)
% LY = str2double(input('Insert Reservoir Width (m) = ','s')); %width of the reservoir (m)
% NX = round(str2double(input('Insert Number of Grid Cells in x direction = ','s')));
% NY = round(str2double(input('Insert Number of Grid Cells in y direction = ','s')));
% PL = str2double(input('Insert Pressure at Left Boundary (Pa) = ','s')); % Pressure at left boundary
% PR = str2double(input('Insert Pressure at Right Boundary (Pa) = ','s')); % Pressure at right boundary 
% PT = str2double(input('Insert Pressure at Top Boundary (Pa) = ','s')); % Pressure at left boundary
% PB = str2double(input('Insert Pressure at Bottom Boundary (Pa) = ','s')); % Pressure at left boundary
% well = str2double(input('Do you want to put wells (Press 1 for Yes, press 0 for No) = ','s'));
% 
% if well == 1
%    nwell = str2double(input('How many wells do you want to put = ','s')); % Number of wells
%    fprintf('You have %g number of grid cell \n',N)
%    fprintf('In which grid cell you want to put wells (from 1 to %g ) = \n', N)
%    for i = 1 : nwell
%     fprintf('Well(%s) : \n',num2str(i))
%     grid(i) = str2double(input('Enter well location = ','s'));
%     if grid(i) > N
%         while grid(i) > N
%         disp('Invalid, well location exceeds number of grid cell')
%         grid(i) = str2double(input('Enter new well location = ','s'));
%         end
%     end
%     pi(i) = str2double(input('Enter productivity index = ','s'));
%     pw(i) = str2double(input('Enter wellbore pressure = ','s'));
%    end
% end
% 
% homogenous = str2double(input('Press 1 for homogenous reservoir, 0 for heterogenous reservoir = ','s')); % case homogenous or heterogenous

%%
NX = 20; %number of grid cells in x direction
NY = 20; %number of grid cell in y directionn
LX = 10; %length of the reservoir (m)
LY = 10; %height of the reservoir (m) 
PL = 10e0; %input boundary condition
PR = 10e0; %input boundary condition
PT = 10e0;
PB = 10e0;

well = 1;
nwell = 2;
grid = [1 NX*NY];
pi = [1000 1000];
pw = [1 0];
homogenous = 1;


%% Initialization
dx = LX/NX; %grid size
x = linspace(dx/2,LX-dx/2,NX); %location of grid center
xi = linspace(0,LX,NX+1); %location of interfaces

dy = LY/NY; %grid size
y = linspace(dy/2,LY-dy/2,NY); %location of grid center
yi = linspace(0,LY,NY+1); %location of interfaces

ceff = 0.1; % effective compressibility in 1/Pa
k = 10^-13; % absolute permeability in m2
por = 0.1;
visw = 10^-3; % water viscosity
densw = 1000; % water density in kg/m3

lambda = zeros(NY,NX);
lambdaharx = zeros(NY,NX+1);
lambdahary = zeros(NY+1,NX);
TX = zeros(NY,NX+1);
TY = zeros(NY+1,NX);

switch homogenous
    case(1)
    lambda = ones(NY,NX); % input lambda
    case(0)
    lambda = rand(NY,NX);
end

[TX,TY,lambdaharx,lambdahary] = computetransmisibility2D(lambda,dx,NX,dy,NY);

%% Pressure Solver

A = zeros(NX*NY,NX*NY);
P = zeros(NX*NY,1);
q = zeros(NX*NY,1);
I = zeros(NY,NX);

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

%% Insert BC
switch well
    case(0) % case no well
    for i = 1 : NX
        for j = 1 : NY
        [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
        if i == 1; % left boundary
        A(I,I) = A(I,I)+ TX(j,i);
        q(I) = q(I)+TX(j,i)*PL;

        elseif i == NX; % right boundary
        A(I,I) = A(I,I)+ TX(j,i+1);
        q(I) = q(I)+TX(j,i+1)*PR;
        
        elseif j == 1
        A(I,I)= A(I,I)+TY(j,i);
        q(I) = q(I)+TY(j,i)*PT;
        
        elseif j == NY
        A(I,I)= A(I,I)+TY(j+1,i);
        q(I) = q(I)+TY(j+1,i)*PB;
        
        end
        end
    end
    
    case(1) % case with wells, constant PI and Pw
    [A,q] = wells2D(pi,pw,lambda,A,q,grid);
end

P = A\q;
P = P';
[UX,UY] = computevelocity2D(lambdaharx,lambdahary,dx,dy,P,NX,NY,i,j);

%% Final Pressure and Velocity Plot 

P = vec2mat(P,NX);
[xx yy] = meshgrid(x,y);
surf(xx,yy,P)
xlabel('x location')
ylabel('y location')
zlabel('Pressure')
title('Incompressible Pressure Solution')

figure
subplot(1,2,1)
surf(UX)
xlabel('x location')
ylabel('y location')
zlabel('Velocity')
xlim([0 NX])
ylim([0 NY])
title('Incompressible Velocity Solution in x Direction')

subplot(1,2,2)
surf(UY)
xlabel('x location')
ylabel('y location')
zlabel('Velocity')
xlim([0 NX])
ylim([0 NY])
title('Incompressible Velocity Solution in y Direction')

end
