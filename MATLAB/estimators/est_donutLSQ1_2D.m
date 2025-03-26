function xest=est_donutLSQ1_2D(photonsi,patternpos,L,fwhm,background)
if nargin<5
    background=0;
end
photonsi=photonsi-background/length(photonsi);
    pi=photonsi/sum(photonsi);
    % eq 2.63
    xest=-1/(1-(L^2*log(2)/fwhm^2))*sum(pi.*patternpos);
end