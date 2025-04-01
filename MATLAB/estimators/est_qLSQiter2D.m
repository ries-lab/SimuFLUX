function xest=est_qLSQiter2D(photonsi,L,probecenter, iter,eps)
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
xest(abs(xest)>L*1.2)=NaN;
end

function dro=it4o(p,xy,L)
x0=xy(1);y0=xy(2);
d13=p(1)-p(3);
d24=p(2)-p(4);

dr(1)=-1/2*((L^2 + 4*(x0^2 + y0^2))*(2*x0*(L + 4*d24*y0) + d13*(L^2 + 4*(x0 - y0)*(x0 + y0))))/(L^3 - 4*L*(x0^2 + y0^2)); 
dr(2)= -1/2*((L^2 + 4*(x0^2 + y0^2))*(2*(L + 4*d13*x0)*y0 + d24*(L^2 - 4*x0^2 + 4*y0^2)))/(L^3 - 4*L*(x0^2 + y0^2));

dro=dr+xy;
end

function dro=it4c(p,xy,L)
x0=xy(1);y0=xy(2);
d13=p(1)-p(3);
d24=p(2)-p(4);
ps=4*p(5)-p(1)-p(2)-p(3)-p(4);

dr(1)=-(((L^2 + 5*(x0^2 + y0^2))*(-(x0*(L^3*(-3 + ps) - 15*d24*L^2*y0 + 5*L*(4 + ps)*(x0^2 + y0^2) + 100*d24*y0*(x0^2 + y0^2))) + ...
     d13*(2*L^4 - 15*L^2*y0^2 + 50*(-x0^4 + y0^4))))/(4*L^5 - 30*L^3*(x0^2 + y0^2) + 100*L*(x0^2 + y0^2)^2)); 
dr(2)= -(((L^2 + 5*(x0^2 + y0^2))*(-(y0*(L^3*(-3 + ps) - 15*d13*L^2*x0 + 5*L*(4 + ps)*(x0^2 + y0^2) + 100*d13*x0*(x0^2 + y0^2))) + ...
     d24*(2*L^4 - 15*L^2*x0^2 + 50*(x0^4 - y0^4))))/(4*L^5 - 30*L^3*(x0^2 + y0^2) + 100*L*(x0^2 + y0^2)^2));

dro=dr+xy;
end

function dro=it6o(p,xy,L)
x0=xy(1);y0=xy(2);
% p=p([2 3 4 5 6 1]);
% p=p([ 1 6 5 4 3 2]);
p1=p(1);p2=p(2);p3=p(3);p4=p(4);p5=p(5);p6=p(6);

dr(1)= -1/4*((4*(2*p1 + p2 - p3 - 2*p4 - p5 + p6)*x0^2 + 4*x0*(L + 2*sqrt(3)*(p2 + p3 - p5 - p6)*y0) + (2*p1 + p2 - p3 - 2*p4 - p5 + p6)*(L^2 - 4*y0^2))*(L^2 + 4*(x0^2 + y0^2)))/(L^3 - 4*L*(x0^2 + y0^2)); 
dr(2)= -1/4*(sqrt(3)*(p2 + p3 - p5 - p6)*(L^4 - 16*x0^4) + 4*(L + 2*(2*p1 + p2 - p3 - 2*p4 - p5 + p6)*x0)*(L^2 + 4*x0^2)*y0 + 8*sqrt(3)*L^2*(p2 + p3 - p5 - p6)*y0^2 + 16*(L + 2*(2*p1 + p2 - p3 - 2*p4 - p5 + p6)*x0)*y0^3 + ... 
    16*sqrt(3)*(p2 + p3 - p5 - p6)*y0^4)/(L^3 - 4*L*(x0^2 + y0^2));
dro=dr+xy;
end

function dro=it6c(p,xy,L)
x0=xy(1);y0=xy(2);

p1=p(1);p2=p(2);p3=p(3);p4=p(4);p5=p(5);p6=p(6);p0=p(7);


dr(1)= -1/12*((3*L^2 + 14*(x0^2 + y0^2))*(-196*(2*p1 + p2 - p3 - 2*p4 - p5 + p6)*x0^4 + 28*x0^3*(L*(-6 - 6*p0 + p1 + p2 + p3 + p4 + p5 + p6) + 14*sqrt(3)*(-p2 - p3 + p5 + p6)*y0) + ...
     2*x0*(3*L^3*(5 - 6*p0 + p1 + p2 + p3 + p4 + p5 + p6) + 35*sqrt(3)*L^2*(p2 + p3 - p5 - p6)*y0 + 14*L*(-6 - 6*p0 + p1 + p2 + p3 + p4 + p5 + p6)*y0^2 + 196*sqrt(3)*(-p2 - p3 + p5 + p6)*y0^3) + ...
     (2*p1 + p2 - p3 - 2*p4 - p5 + p6)*(9*L^4 - 70*L^2*y0^2 + 196*y0^4)))/(9*L^5 - 70*L^3*(x0^2 + y0^2) + 196*L*(x0^2 + y0^2)^2); 
 dr(2)= -1/12*((3*L^2 + 14*(x0^2 + y0^2))*(9*sqrt(3)*L^4*(p2 + p3 - p5 - p6) + 6*L^3*(5 - 6*p0 + p1 + p2 + p3 + p4 + p5 + p6)*y0 + 70*L^2*x0*(sqrt(3)*(-p2 - p3 + p5 + p6)*x0 + (2*p1 + p2 - p3 - 2*p4 - p5 + p6)*y0) + ...
     28*L*(-6 - 6*p0 + p1 + p2 + p3 + p4 + p5 + p6)*y0*(x0^2 + y0^2) + 196*(sqrt(3)*(p2 + p3 - p5 - p6)*x0^4 + 2*(-2*p1 - p2 + p3 + 2*p4 + p5 - p6)*x0^3*y0 + 2*(-2*p1 - p2 + p3 + 2*p4 + p5 - p6)*x0*y0^3 + ...
       sqrt(3)*(-p2 - p3 + p5 + p6)*y0^4)))/(9*L^5 - 70*L^3*(x0^2 + y0^2) + 196*L*(x0^2 + y0^2)^2);

dro=dr+xy;
end

function dro=it3c(p,xy,L)
x0=xy(1);y0=xy(2);
% p=p([ 1 6 5 4 3 2]);
p1=p(1);p2=p(2);p3=p(3);p0=p(4);

dr(1)= -1/12*((3*L^2 + 16*(x0^2 + y0^2))*(9*L^4*(2*p1 - p2 - p3) + 12*L^3*(2 - 3*p0 + p1 + p2 + p3)*x0 + 64*L^2*y0*(sqrt(3)*(p2 - p3)*x0 + (-2*p1 + p2 + p3)*y0) - 64*L*(3 + 3*p0 - p1 - p2 - p3)*x0*(x0^2 + y0^2) - ...
     256*(x0^2 + y0^2)*((2*p1 - p2 - p3)*x0^2 + 2*sqrt(3)*(p2 - p3)*x0*y0 + (-2*p1 + p2 + p3)*y0^2)))/(9*L^5 - 64*L^3*(x0^2 + y0^2) + 256*L*(x0^2 + y0^2)^2); 
dr(2)= -1/12*((3*L^2 + 16*(x0^2 + y0^2))*(9*sqrt(3)*L^4*(p2 - p3) + 12*L^3*(2 - 3*p0 + p1 + p2 + p3)*y0 + 64*L*(-3 - 3*p0 + p1 + p2 + p3)*y0*(x0^2 + y0^2) + ...
     256*(sqrt(3)*(p2 - p3)*x0^4 + 2*(-2*p1 + p2 + p3)*x0^3*y0 + 2*(-2*p1 + p2 + p3)*x0*y0^3 + sqrt(3)*(-p2 + p3)*y0^4) - 64*L^2*x0*(-2*p1*y0 + p3*(-(sqrt(3)*x0) + y0) + p2*(sqrt(3)*x0 + y0))))/...
   (9*L^5 - 64*L^3*(x0^2 + y0^2) + 256*L*(x0^2 + y0^2)^2);
dro=dr+xy;
end

function dro=it3o(p,xy,L)
x0=xy(1);y0=xy(2);
p1=p(1);p2=p(2);p3=p(3);

dr(1)= -1/4*((4*(2*p1 - p2 - p3)*x0^2 + 4*x0*(L + 2*sqrt(3)*(p2 - p3)*y0) + (2*p1 - p2 - p3)*(L^2 - 4*y0^2))*(L^2 + 4*(x0^2 + y0^2)))/(L^3 - 4*L*(x0^2 + y0^2));
dr(2)= -1/4*((L^2 + 4*(x0^2 + y0^2))*(sqrt(3)*L^2*(p2 - p3) + 4*L*y0 - 4*(sqrt(3)*(p2 - p3)*x0^2 + 2*(-2*p1 + p2 + p3)*x0*y0 + sqrt(3)*(-p2 + p3)*y0^2)))/(L^3 - 4*L*(x0^2 + y0^2));
dro=dr+xy;
end