function [q,dqdP] = computewellflux1D(pw,lambda,pi,grid,P)
N = length(P);
q = zeros(N,1);
dqdP = zeros(N,1);
nwell= length(pw);
    for i = 1 : nwell
        w = grid(i);
        q(w)=pi(i)*lambda(w)*(pw(i)-P(w));
        dqdP(w)= -pi(i)*lambda(w);
    end

end