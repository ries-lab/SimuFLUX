function xest=est_zline(photonsi,L,iter,eps)
if nargin<4
    eps=.1; %nm
end
if nargin<3
    iter=15;
end
phot1=photonsi(1:2);
phot2=photonsi(3:4);
if length(photonsi)==5
   phot1(end+1)=photonsi(5);
   phot2(end+1)=photonsi(5);
end

c=[2.1415   -3.2122    2.0062   -0.4677]*1e3;
fr=phot2(2)/(phot2(1)+phot2(2));
xest=c(1)*fr^3+c(2)*fr^2+c(3)*fr+c(4);
if abs(xest)<L(1)/2*0.75
% xest1=est_qLSQiter1D(phot2,L(2),iter,eps,0);
    xest=est_qLSQiter1D(phot1,L(1),iter,eps,xest);
end
end
