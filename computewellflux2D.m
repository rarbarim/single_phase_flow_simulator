function q = computewellflux2D(pw,lambda,pi,grid,P)
N = length(P);
q = zeros(N,1);
nwell= length(pw);
    for w = 1 : nwell
        q(grid(w))=pi(w)*lambda(grid(w))*(pw(w)-P(grid(w)));
    end

end