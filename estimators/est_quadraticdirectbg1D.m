function xest=est_quadraticdirectbg1D(photonsi,L)
if length(photonsi)>3
    i1=3;i2=1;i0=5;
else
    i1=1;i2=2;i0=3;
end
p=photonsi([i1 i2 i0]);%for testing:
xest=(-(L*p(1)) + L*p(2))/(4*(2*p(3) - p(1) - p(2)));
xest(abs(xest)>L)=NaN;
end
