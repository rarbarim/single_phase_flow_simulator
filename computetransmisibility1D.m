function [T,lambdahar] = computetransmisibility1D(lambda,dx,N)

lambdahar(1) = lambda(1);
lambdahar(2:N) = 2*lambda(2:N).*lambda(1:N-1)./(lambda(2:N)+lambda(1:N-1));
lambdahar(N+1) = lambda(N);

T(1)=lambdahar(1)./((dx^2)/2);
T(2:N)=lambdahar(2:N)./(dx^2);
T(N+1)=lambdahar(N+1)./((dx^2)/2);

end