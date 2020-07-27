function [P,UX,UY] = slightly_compressible_flow_solver2D();
% 2D eliptic FD and FV solver
% slightly compresisble flow
% |---*---|---*---|---*---|
% Note : number of interface = number of cells + 1

% close all
% clear all

%% User input  

% disp('Welcome in 1D Slightly Compressible Flow Simulator by RB Arbarim - 4573900')
% disp(' ')
% LX = str2double(input('Insert Reservoir Length (m) = ','s')); %length of the reservoir (m)
% LY = str2double(input('Insert Reservoir Width (m) = ','s')); %width of the reservoir (m)
% NX = round(str2double(input('Insert Number of Grid Cells in x direction = ','s')));
% NY = round(str2double(input('Insert Number of Grid Cells in y direction = ','s')));
% PL = str2double(input('Insert Pressure at Left Boundary (Pa) = ','s')); % Pressure at left boundary
% PR = str2double(input('Insert Pressure at Right Boundary (Pa) = ','s')); % Pressure at right boundary 
% PT = str2double(input('Insert Pressure at Top Boundary (Pa) = ','s')); % Pressure at left boundary
% PB = str2double(input('Insert Pressure at Bottom Boundary (Pa) = ','s')); % Pressure at left boundary
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

LX = 10;
LY = 10; % LY becomes 1 in 1D
NX = 20;
NY = 20; % NY becomes 1 in 1D
PL = 10e0;;
PR = 10e0;;
PT = 10e0;;
PB = 10e0;;
t = 10e5;
nt = 10;
P0 = 1;
cf = 1;
cr = 1;
well = 1;
nwell = 2;
grid = [1 NX*NY];
pi = [1000 1000];
pw = [1 0];
homogenous = 1;
solution = 0;

%% Initialization
dx = LX/NX; %grid size
x = linspace(dx/2,LX-dx/2,NX); %location of grid center
xi = linspace(0,LX,NX+1); %location of interfaces

dy = LY/NY; %grid size
y = linspace(dy/2,LY-dy/2,NY); %location of grid center
yi = linspace(0,LY,NY+1); %location of interfaces

dt = t/nt; % time grid size 
ceff = 10e-6; % effective compressibility in 1/Pa
k = 10^-13; % absolute permeability in m2
por = 0.1;
visw = 10^-3; % water viscosity
densw = 1000; % water density in kg/m3
[xx yy] = meshgrid(x,y);

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

%% Time Loop
P = zeros(NX*NY,1);

for Ndt = 1 : round(nt)
A = zeros(NX*NY,NX*NY);
q = zeros(NX*NY,1);
I = zeros(NY,NX);

for i = 1 : NX
    for j = 1:NY
            gridindex(j,i) = (j-1)*NX+i;
            I = (j-1)*NX+i;
            IE = (j-1)*NX+i+1;
            IW = (j-1)*NX+i-1;
            IS = (j)*NX+i;
            IN = (j-2)*NX+i;
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
%     end
end

C = eye(NX*NY)*(ceff*por/dt);

% solution=0; %0 for explicit, 1 for implicite
switch solution
    case(1)
    P = C\(q-A*P+C*P);
    case(0)
    P = (C+A)\(q+C*P);
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
title('Slightly Compressible Pressure Solution')
hold off

%% Plot solution

P = vec2mat(P,NX);
figure
surf(xx,yy,P)
xlabel('x location')
ylabel('y location')
zlabel('Pressure')
title('Slightly Compressible Pressure Solution')

figure
subplot(1,2,1)
surf(UX)
xlabel('x location')
ylabel('y location')
zlabel('Velocity')
xlim([0 NX])
ylim([0 NY])
title('Slightly Compressible Velocity Solution in x Direction')

subplot(1,2,2)
surf(UY)
xlabel('x location')
ylabel('y location')
zlabel('Velocity')
xlim([0 NX])
ylim([0 NY])
title('Slightly Compressible Velocity Solution in y Direction')

end