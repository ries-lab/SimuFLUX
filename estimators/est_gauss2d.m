function xest=est_gauss2d(photonsi,patternpos,L,sigma,iscenter)
    numrim=size(patternpos,1)-iscenter;
    pi=photonsi/sum(photonsi);
    Ls2=L^2/8/sigma^2;
    % eq 2.63
    xest=(numrim+exp(Ls2))/Ls2/numrim*sum(pi.*patternpos);
end