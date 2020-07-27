% this is the main file which forward the user into desired 1D flow system 

clear all
close all
format short

%% forward
disp('Welcome in 1D Flow Simulator by RB Arbarim - 4573900')
disp(' ')
disp('Press 1 for Incompressible Flow Solver')
disp('Press 2 for Slightly Compressible Flow Solver')
disp('Press 3 for Compressible Flow Solver')

% system = str2double(input('Choose Option = ','s')); 
system = 1; % 1 for incompressible, 2 for slightly compressible, 3 for compressible
name = cell(3,1);
name = {'Incompressible Flow Solver' 'Slightly Compressible Flow Solver' 'Compressible Flow Solver'};
fprintf('Please wait,you will be redirected to %s \n', cell2mat(name(system)))

switch system
    case(1)
        [P,U] = incompressible_flow_solver1D();
    case(2)
        [P,U] = slightly_compressible_flow_solver1D();
        [P2] = slightly_compressible_flow_solver1D_analytical(); % validation 
    case(3)
        [P,U] = compressible_flow_solver1D();
end

disp('')
disp('VOILA!! see the figure, thanks for using this simulator')