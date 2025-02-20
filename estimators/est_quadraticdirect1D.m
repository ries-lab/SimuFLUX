function xest=est_quadraticdirect1D(photonsi,L)
if length(photonsi)>3
    i1=3;i2=1;i0=5;
else
    i1=1;i2=2;i0=3;
end
maxrange=150;
    pi=photonsi/(photonsi(i1)+photonsi(i2));
    % eq 2.63
    xest1=(L-sqrt(L^2*(1-(pi(i1)-pi(i2))^2)))/(2*(pi(i1)-pi(i2)));
    xest2=(L+sqrt(L^2*(1-(pi(i1)-pi(i2))^2)))/(2*(pi(i1)-pi(i2)));
    if xest2>maxrange || sum((pnorm(L,xest1)-pi([i1,i2,i0])).^2)< sum((pnorm(L,xest2)-pi([i1,i2,i0])).^2)
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