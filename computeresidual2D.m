function residual = computeresidual2D(rhon,phin,rho,phi,TX,TY,dt,P,q,NX,NY);

A = zeros(NX*NY,NX*NY);
%         q = zeros(N,1);

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
      
%     for i = 1:N
%         if i > 1 % there is left neighbor
%             A(i,i)= T(i);
%             A(i,i-1)= -T(i);
%         end
%         if i < N % there is right neighbor
%             A(i,i)= A(i,i)+T(i+1);
%             A(i,i+1)= -T(i+1);
%         end
%     end
    
residual = (rho.*q)-(((rho.*phi)-(rhon.*phin))/dt)-A*P;

 
end 