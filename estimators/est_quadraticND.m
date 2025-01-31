function xest=est_quadraticND(photonsi,patternpos)
    pi=photonsi/sum(photonsi);
    % eq 2.63
    xest=-sum(pi'.*patternpos);
end