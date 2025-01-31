function xest=est_donut2d(photonsi,patternpos,L,fwhm,iscenter)
    pi=photonsi/sum(photonsi);
    % eq 2.63
    xest=-1/(1-(L^2*log(2)/fwhm^2))*sum(pi.*patternpos);
end