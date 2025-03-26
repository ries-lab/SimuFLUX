function xest=est_qMLE1D(photonsi,L)
if length(photonsi)>3
    i1=3;i2=1;i0=5;
else
    i1=1;i2=2;i0=3;
end
    maxrange=L;
    pi=photonsi/(photonsi(i1)+photonsi(i2));
    % eq 2.63
    xest1=-L/2+L/(1+sqrt(photonsi(i2)/photonsi(i1)));
    xest2=-L/2+L/(1-sqrt(photonsi(i2)/photonsi(i1)));
    
    p1=pnorm(L,xest1);p2=pnorm(L,xest2);
    LL1=loglikelihood(p1/sum(p1),photonsi([i1 i2 i0]));
    LL2=loglikelihood(p2/sum(p2),photonsi([i1 i2 i0]));

    if LL1>LL2
        xest=xest1;
    else
        xest=xest2;
    end
    % xest=min(xest,150);xest=max(xest,-150); %deviation of Gaussian from Donut
    xest(abs(xest)>maxrange)=NaN; %deviation of Gaussian from Donut
end

function po=pnorm(L,x0)
    po=[1/2+(2*L*x0)/(L^2+4*x0^2);1/2-(2*L*x0)/(L^2+4*x0^2);2*x0^2/(L^2+4*x0^2)];
end

function LL=loglikelihood(pi,ni)
LL=sum(ni.*log(pi));
end