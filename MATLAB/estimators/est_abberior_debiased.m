function xest=est_abberior_debiased(photonsi,patternpos,coefficients,patGeoFactor)
    nmtorel=360*patGeoFactor;
    % nmtorel=1;
    pi=photonsi/sum(photonsi);
    xestd=-sum(pi.*patternpos);
    ch=coefficients(1,:,1);
    xa=sqrt(sum(xestd.^2))/nmtorel;
    cf=1+xa*ch(3)+xa^2*ch(2)+xa^3*ch(1)
    xest=cf*xestd;
    asdf
end