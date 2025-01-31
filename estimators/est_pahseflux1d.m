function xest=est_pahseflux1d(photonsi,coord,L)
xest=[0 0 0];
    % ph=photonsi/sum(photonsi);
    % eq 2.63
    xest(coord)=L/(1+sqrt(photonsi(end)/photonsi(1)))-L/2;
end