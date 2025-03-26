function xest=est_qLSQ1_2D(photonsi,patternpos)
    pi=photonsi/sum(photonsi);
    % eq 2.63
    xest=-sum(pi.*patternpos);
end