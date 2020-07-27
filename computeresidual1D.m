% function R = computeresidual1D(rhon,phin,rho,phi,T,dt,P,q,N);
function [R,J] = computeresidual1D(rhon,phin,rho,drho,phi,dphi,T,dt,P,q,lambda,pi,N,grid);

A = zeros(N,N);
%         q = zeros(N,1);

    for i = 1:N
        if i > 1 % there is left neighbor
            A(i,i)= T(i);
            A(i,i-1)= -T(i);
        end
        if i < N % there is right neighbor
            A(i,i)= A(i,i)+T(i+1);
            A(i,i+1)= -T(i+1);
        end
    end
    
R = -((rho.*q)-(((rho.*phi)-(rhon.*phin))/dt)-A*P);

%% Jacobian
    vec = (dphi.*rho + drho.*phi)/dt;
    C = diag(vec);

    W = zeros(N,N); % only non zero at well location
    
   for i = 1:length(pi);
       W(grid(i),grid(i))= W(grid(i),grid(i))+lambda(grid(i))*pi(i)*rho(grid(i)); % from class
   end
    
    J = A+C+W; detJ = det(J);
 
end 