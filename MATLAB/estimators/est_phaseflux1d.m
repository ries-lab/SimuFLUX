function xest=est_phaseflux1d(photonsi,L)
% xest=[0 0 0];
    % ph=photonsi/sum(photonsi);
    % eq 2.63
    xest=L/(1+sqrt(photonsi(end)/photonsi(1)))-L/2;
end