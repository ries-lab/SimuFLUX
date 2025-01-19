function xest=estimators(estimatorname, photons, patternpos, varargin)
switch estimatorname
    case 'quadraticND'
        xest=positionestimatequad(photons,patternpos);
    case 'donut2D'
        xest=positionestimatedonut(photons,patternpos,varargin{:});
    case 'phaseflux1D'
        xest=positionestimate1D(photons,patternpos,varargin{:});
    case 'donut1D'
        xest=positionestimatedonut1D(photons,patternpos,varargin{:});
    case 'gauss2D'
        xest=positionestimategauss2D(photons,patternpos,varargin{:});
    otherwise
        disp('estimator not found')
end

end



function xest=positionestimatequad(photonsi,patternpos)
    pi=photonsi/sum(photonsi);
    % eq 2.63
    xest=-sum(pi'.*patternpos);
end

function xest=positionestimatedonut(photonsi,patternpos,L,fwhm,probecenter)
    pi=photonsi/sum(photonsi);
    % eq 2.63
    xest=-1/(1-(L^2*log(2)/fwhm^2))*sum(pi.*patternpos);
end

function xest=positionestimate1D(photonsi,coord,L)
xest=[0 0 0];
    % ph=photonsi/sum(photonsi);
    % eq 2.63
    xest(coord)=L/(1+sqrt(photonsi(end)/photonsi(1)))-L/2;
end

function xest=positionestimatedonut1D(photonsi,patternpos,coord,sigma)
xest=[0 0 0];

ph=photonsi([1 2 3]); pp=patternpos([1 2 3]);
 pis=ph/sum(ph);
 L=2*patternpos(end);
    % ph=photonsi/sum(photonsi);
    % eq 2.63
    xest(coord)=-1/(1-(L^2/sigma^2))*sum(pis'.*pp)*2/4;
    xest(coord)=-sum(pis'.*pp)/2;
end



function xest=positionestimategauss2D(photonsi,patternpos,L,sigma,iscenter)
    numrim=size(patternpos,1)-iscenter;
    pi=photonsi/sum(photonsi);
    Ls2=L^2/8/sigma^2;
    % eq 2.63
    xest=(numrim+exp(Ls2))/Ls2/numrim*sum(pi.*patternpos);
end