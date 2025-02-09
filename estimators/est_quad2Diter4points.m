function xest=est_quad2Diter4points(photonsi,L,iter,eps)
%Eilers 2.64, k=1
if nargin<4
    eps=.1; %nm
end
if nargin<3
    iter=15;
end
if length(photonsi)==4
    itfun=@iteration;
else
    itfun=@iterationcenter;
end
    pi=photonsi/sum(photonsi);
    xest=[0,0];
    for k=1:iter
        xo=xest;
        xest=itfun(pi,xest,L);
        if sum((xest-xo).^2)<eps^2
            break
        end
    end    
end

function dro=iteration(p,xy,L)
x0=xy(1);y0=xy(2);
d13=p(2)-p(4);
d24=p(3)-p(1);

dr(1)= (((L^2 + 4 *(x0^2 + y0^2)) *(-2 * x0 * (L + 4 * d13 * y0) + ...
    d24 * (L^2 + 4 * (x0 - y0)* (x0 + y0))))/( 2 *(L^3 - 4 *L * (x0^2 + y0^2))));
dr(2)=-(((L^2 + 4 * (x0^2 + y0^2)) * (2 * (L - 4 * d24 * x0) * y0 + ...
     d13 *(L^2 - 4 * x0^2 + 4 *y0^2)))/(2 * (L^3 - 4 *L *(x0^2 + y0^2))));
dro=dr+xy;
end

function dro=iterationcenter(p,xy,L)
x0=xy(1);y0=xy(2);
d13=p(2)-p(4);
d24=p(3)-p(1);
ps=4*p(5)-p(1)-p(2)-p(3)-p(4);

dr(1)=((L^2 + 5 *(x0^2 + y0^2)) * (x0 * (L^3 *(-3 + ps) - 15 *d13 *L^2 * y0 + ...
         5 *L *(4 + ps) *(x0^2 + y0^2) + 100 * d13 * y0 * (x0^2 + y0^2)) + ... 
      d24 * (2 * L^4 - 15  * L^2 * y0^2 + 50 * (-x0^4 + y0^4))))/(4 * L^5 - ...
    30 * L^3 *(x0^2 + y0^2) + 100 * L * (x0^2 + y0^2)^2);
       
       
dr(2)= -(((L^2 +  5 *(x0^2 + y0^2)) *(-y0 * (L^3 *(-3 + ps) + 15* d24 *L^2 *x0 + ...
           5 *L *(4 + ps)* (x0^2 + y0^2) - 100 *d24 *x0 *(x0^2 + y0^2)) + ...
        d13* (2 *L^4 - 15 *L^2* x0^2 + 50 *(x0^4 - y0^4))))/(4 *L^5 - ...
      30* L^3 *(x0^2 + y0^2) + 100 * L * (x0^2 + y0^2)^2));

dro=dr+xy;
end