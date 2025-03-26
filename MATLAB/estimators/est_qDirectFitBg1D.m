function xest=est_qDirectFitBg1D(photonsi,L,probecenter)
if length(photonsi)>3
    i1=3;i2=1;i0=5;
elseif length(photonsi)==3
    i1=1;i2=2;i0=3;
elseif length(photonsi)==2
    xest=q2point(photons,L);
    return
end
p=photonsi([i1 i2 i0]);%for testing:
% xest=(-(L*p(1)) + L*p(2))/(4*(2*p(3) - p(1) - p(2)));
xest=(L*(-p(1) + p(2)))/(8*p(3) - 4*(p(1) + p(2)));
relmax=1;
xest(abs(xest)>relmax*L)=relmax*L; %extrapolation does not work well
end

function xest=q2point(photons,L)
maxrange=150;
xest1=(L*(sqrt(p1) - sqrt(p2)))/(2*(sqrt(p1) + sqrt(p2)));
xest2=(L*(sqrt(p1) + sqrt(p2))^2)/(2*(p1 - p2));
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