function xest=est_GaussLSQ1_2D(photonsi,patternpos,L,sigma,iscenter)
    pi=photonsi/sum(photonsi);
    Ls2=L^2/8/sigma^2;
    % eq 2.63
    if iscenter
        numrim=size(patternpos,1)-iscenter;
        xest=(numrim+exp(Ls2))/Ls2/numrim*sum(pi.*patternpos);
    else
        xest=sum(pi.*patternpos)/Ls2;
    end
end