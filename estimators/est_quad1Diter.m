function xest=est_quad1Diter(photonsi,L,iter,eps)
%Eilers 2.64, k=1
if nargin<4
    eps=.1; %nm
end
if nargin<3
    iter=15;
end
if length(photonsi)==2
    itfun=@iteration;
else
    itfun=@iterationcenter;
end
pi=photonsi/sum(photonsi);
xest=0;
for k=1:iter
    xo=xest;
    xest=itfun(pi,xest,L);
    if sum((xest-xo).^2)<eps^2
        break
    end
end 
% xest=max(min(xest,L),-L); %avoid crazy high numbers
xest(abs(xest)>1.3*L)=NaN;
end

function dro=iteration(p,x0,L)
d12=p(1)-p(2);
dr=((L^2 + 4 *x0^2)* (-4 *L *x0 + d12* (L^2 + 4 *x0^2)))/(4 *(L^3 - 4 *L *x0^2));
dro=dr+x0;
end

function dro=iterationcenter(p,x0,L)
%p(1) at -L/2, p(2) at L/2, p(3) at 0
d12=p(1)-p(2);
pt=p(1)+p(2)-2*p(3);
dr=((L^2 + 6 *x0^2) *(-L^3 *(3 + pt) *x0 - 6 *L *(-4 + pt) *x0^3 + ...
   d12 *(L^4 - 36 *x0^4)))/(4 *(L^5 - 9 *L^3 *x0^2 + 36 *L *x0^4));
dro=dr+x0;
end