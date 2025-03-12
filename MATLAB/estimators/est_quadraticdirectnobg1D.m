function xest=est_quadraticdirectnobg1D(photonsi,L,probecenter)
if length(photonsi)>3
    p1=photonsi(3);p2=photonsi(1);
else
    p1=photonsi(2);p2=photonsi(1);
end
if ~probecenter
    xest=(L*(sqrt(p1) - sqrt(p2)))/(2*(sqrt(p1) + sqrt(p2)));
    return
end
p0=photonsi(end);
maxrange=150;
xest1=(L*(sqrt(p1) - sqrt(p2)))/(2*(sqrt(p1) + sqrt(p2)));
xest2=(L*(sqrt(p1) + sqrt(p2))^2)/(2*(p1 - p2));
    if xest2>maxrange || sum((pnorm(L,xest1)-[p1;p2;p0]/(p1+p2)).^2)<= sum((pnorm(L,xest2)-[p1;p2;p0]/(p1+p2)).^2)
        xest=xest1;
    else
        xest=xest2;
    end
    % xest=min(xest,150);xest=max(xest,-150); %deviation of Gaussian from Donut
    relmax=0.75;
    xest(abs(xest)>relmax*L)=relmax*L; %deviation of Gaussian from Donut
end


function po=pnorm(L,x0)
    po=[1/2+(2*L*x0)/(L^2+4*x0^2);1/2-(2*L*x0)/(L^2+4*x0^2);2*x0^2/(L^2+4*x0^2)];
end