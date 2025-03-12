function xest=est_quad2Diter(photonsi,L,probecenter, iter,eps)
%Eilers 2.64, k=1
if nargin<5
    eps=.1; %nm
end
if nargin<4
    iter=15;
end
if nargin<3 || isempty(probecenter)
    probecenter=mod(length(photonsi),2);
end
orbitpoints=length(photonsi)-probecenter;
if probecenter
    ct="c";
else
    ct="o";
end
itfun=str2func("it"+orbitpoints+ct);
pi=photonsi/sum(photonsi);
xest=[0,0];
for k=1:iter
    xo=xest;
    xest=itfun(pi,xest,L);
    if sum((xest-xo).^2)<eps^2
        break
    end
end    
% xest=max(min(xest,L),-L); %avoid crazy high numbers
xest(abs(xest)>L)=NaN;
end

function dro=it4o(p,xy,L)
x0=xy(1);y0=xy(2);
d13=p(2)-p(4);
d24=p(3)-p(1);

dr(1)= (((L^2 + 4 *(x0^2 + y0^2)) *(-2 * x0 * (L + 4 * d13 * y0) + ...
    d24 * (L^2 + 4 * (x0 - y0)* (x0 + y0))))/( 2 *(L^3 - 4 *L * (x0^2 + y0^2))));
dr(2)=-(((L^2 + 4 * (x0^2 + y0^2)) * (2 * (L - 4 * d24 * x0) * y0 + ...
     d13 *(L^2 - 4 * x0^2 + 4 *y0^2)))/(2 * (L^3 - 4 *L *(x0^2 + y0^2))));
dro=dr+xy;
end

function dro=it4c(p,xy,L)
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

function dro=it6o(p,xy,L)
x0=xy(1);y0=xy(2);
% p=p([ 1 6 5 4 3 2]);
p1x=p(1);p2x=p(2);p3x=p(3);p4x=p(4);p5x=p(5);p6x=p(6);

dr(1)= -1/4*((L^2 + 4*(x0^2 + y0^2))*(L^2*(2*p1x + p2x - p3x - 2*p4x - p5x + p6x) + 4*L*x0 + ...
      4*(-(p3x*x0^2) - 2*p4x*x0^2 - p5x*x0^2 + p6x*x0^2 + 2*sqrt(3)*p3x*x0*y0 - 2*sqrt(3)*p5x*x0*y0 - ...
        2*sqrt(3)*p6x*x0*y0 + p3x*y0^2 + 2*p4x*y0^2 + p5x*y0^2 - p6x*y0^2 + 2*p1x*(x0^2 - y0^2) + ...
        p2x*(x0^2 + 2*sqrt(3)*x0*y0 - y0^2))))/(L^3 - 4*L*(x0^2 + y0^2));
dr(2)=-1/4*((L^2 + 4*(x0^2 + y0^2))*(sqrt(3)*L^2*(p2x + p3x - p5x - p6x) + 4*L*y0 + ...
      4*(sqrt(3)*p5x*x0^2 + sqrt(3)*p6x*x0^2 + 4*p1x*x0*y0 - 4*p4x*x0*y0 - 2*p5x*x0*y0 + 2*p6x*x0*y0 - sqrt(3)*p5x*y0^2 - ...
        sqrt(3)*p6x*y0^2 + p3x*(-(sqrt(3)*x0^2) - 2*x0*y0 + sqrt(3)*y0^2) + ...
        p2x*(-(sqrt(3)*x0^2) + 2*x0*y0 + sqrt(3)*y0^2))))/(L^3 - 4*L*(x0^2 + y0^2));
dro=dr+xy;
end

function dro=it6c(p,xy,L)
x0=xy(1);y0=xy(2);
p1x=p(1);p2x=p(2);p3x=p(3);p4x=p(4);p5x=p(5);p6x=p(6);p7x=p(7);

dr(1)= -1/12*((3*L^2 + 14*(x0^2 + y0^2))*(9*L^4*(2*p1x + p2x - p3x - 2*p4x - p5x + p6x) + ...
      6*L^3*(5 + p1x + p2x + p3x + p4x + p5x + p6x - 6*p7x)*x0 + 70*L^2*y0*(sqrt(3)*(p2x + p3x - p5x - p6x)*x0 + ...
        (-2*p1x - p2x + p3x + 2*p4x + p5x - p6x)*y0) + 28*L*(-6 + p1x + p2x + p3x + p4x + p5x + p6x - 6*p7x)*x0*...
       (x0^2 + y0^2) + 196*((-2*p1x - p2x + p3x + 2*p4x + p5x - p6x)*x0^4 + 2*sqrt(3)*(-p2x - p3x + p5x + p6x)*x0^3*y0 + ...
        2*sqrt(3)*(-p2x - p3x + p5x + p6x)*x0*y0^3 + (2*p1x + p2x - p3x - 2*p4x - p5x + p6x)*y0^4)))/...
    (9*L^5 - 70*L^3*(x0^2 + y0^2) + 196*L*(x0^2 + y0^2)^2);

dr(2)= -1/12*((3*L^2 + 14*(x0^2 + y0^2))*(9*sqrt(3)*L^4*(p2x + p3x - p5x - p6x) + ...
      6*L^3*(5 + p1x + p2x + p3x + p4x + p5x + p6x - 6*p7x)*y0 + 70*L^2*x0*(sqrt(3)*(-p2x - p3x + p5x + p6x)*x0 + ...
        (2*p1x + p2x - p3x - 2*p4x - p5x + p6x)*y0) + 28*L*(-6 + p1x + p2x + p3x + p4x + p5x + p6x - 6*p7x)*y0*...
       (x0^2 + y0^2) + 196*(sqrt(3)*(p2x + p3x - p5x - p6x)*x0^4 + 2*(-2*p1x - p2x + p3x + 2*p4x + p5x - p6x)*x0^3*y0 +... 
        2*(-2*p1x - p2x + p3x + 2*p4x + p5x - p6x)*x0*y0^3 + sqrt(3)*(-p2x - p3x + p5x + p6x)*y0^4)))/...
    (9*L^5 - 70*L^3*(x0^2 + y0^2) + 196*L*(x0^2 + y0^2)^2);
dro=dr+xy;
end

function dro=it3c(p,xy,L)
x0=xy(1);y0=xy(2);
% p=p([ 1 6 5 4 3 2]);
p1x=p(1);p2x=p(2);p3x=p(3);p4x=p(4);
dr(1)=-1/12*((3*L^2 + 16*(x0^2 + y0^2))*(9*L^4*(2*p1x - p2x - p3x) + 12*L^3*(2 + p1x + p2x + p3x - 3*p4x)*x0 + ...
      64*L^2*y0*(sqrt(3)*(p2x - p3x)*x0 + (-2*p1x + p2x + p3x)*y0) + 64*L*(-3 + p1x + p2x + p3x - 3*p4x)*x0*...
       (x0^2 + y0^2) - 256*(x0^2 + y0^2)*((2*p1x - p2x - p3x)*x0^2 + 2*sqrt(3)*(p2x - p3x)*x0*y0 + ...
        (-2*p1x + p2x + p3x)*y0^2)))/(9*L^5 - 64*L^3*(x0^2 + y0^2) + 256*L*(x0^2 + y0^2)^2);
 dr(2)=-1/12*((3*L^2 + 16*(x0^2 + y0^2))*(9*sqrt(3)*L^4*(p2x - p3x) + 12*L^3*(2 + p1x + p2x + p3x - 3*p4x)*y0 + ...
      64*L^2*x0*(sqrt(3)*(-p2x + p3x)*x0 - (-2*p1x + p2x + p3x)*y0) + 64*L*(-3 + p1x + p2x + p3x - 3*p4x)*y0*...
       (x0^2 + y0^2) + 256*(sqrt(3)*(p2x - p3x)*x0^4 + 2*(-2*p1x + p2x + p3x)*x0^3*y0 + 2*(-2*p1x + p2x + p3x)*x0*y0^3 + ...
        sqrt(3)*(-p2x + p3x)*y0^4)))/(9*L^5 - 64*L^3*(x0^2 + y0^2) + 256*L*(x0^2 + y0^2)^2);
dro=dr+xy;
end

function dro=it3o(p,xy,L)
x0=xy(1);y0=xy(2);
% p=p([ 1 6 5 4 3 2]);
p1x=p(1);p2x=p(2);p3x=p(3);p4x=p(4);

dr(1)=-1/4*((4*(2*p1x - p2x - p3x)*x0^2 + 4*x0*(L + 2*sqrt(3)*(p2x - p3x)*y0) + (2*p1x - p2x - p3x)*(L^2 - 4*y0^2))*...
     (L^2 + 4*(x0^2 + y0^2)))/(L^3 - 4*L*(x0^2 + y0^2));
dr(2)=-1/4*((L^2 + 4*(x0^2 + y0^2))*(sqrt(3)*L^2*(p2x - p3x) + 4*L*y0 - ...
      4*(sqrt(3)*(p2x - p3x)*x0^2 + 2*(-2*p1x + p2x + p3x)*x0*y0 + sqrt(3)*(-p2x + p3x)*y0^2)))/...
    (L^3 - 4*L*(x0^2 + y0^2));
dro=dr+xy;
end