
function xest=est_donut1d(photonsi,patternpos,L,sigma,iscenter)
% xest=[0 0 0];

ph=photonsi;
 pis=ph/sum(ph);
 % L=2*patternpos(end);
    % ph=photonsi/sum(photonsi);
    % eq 2.63
    xest=-1/(1-(L^2/sigma^2))*sum(pis.*patternpos)*2/4;
    % xest=-sum(pis'.*pp)/2;
end