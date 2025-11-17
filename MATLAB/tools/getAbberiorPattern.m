function [parg,esth]=getAbberiorPattern(itr,seq)
if contains(itr.Mode.epsf,'focFldRing')
    phasemask="tophat";
    sigma_est_ph=130;
elseif contains(itr.Mode.epsf,'focFldVortex')
    phasemask="vortex";
    sigma_est_ph = 190; % Default value for sigma estimation
end
if contains(itr.Mode.modulated,'phl')&& contains(itr.Mode.pattern,'hexagon') 
    pinholeorbit=true;
else
    pinholeorbit=false;
end

% scan patterns
if itr.ccrLimit==-1
    probecenter=false;
else
    probecenter=true;
end
L=itr.patGeoFactor*360; %nm
arg2D={'makepattern','orbitscan'};
dim=[1,2];
switch itr.Mode.pattern
    case 'hexagon'
        arg=arg2D;
        patternpoints=6;      
    case 'square'
        arg=arg2D;    
        patternpoints=4;   
    case 'triangle'
        arg=arg2D; 
        patternpoints=3;   
    case 'zline'
        dim=3;
        patternpoints=length(L)*2;
        patternpos=zeros(patternpoints,3); patternpos(:,3)=[-L(1), L(1), -L(2), L(2)]/2;
        if probecenter %probecenter, pattern points argument ignored if not makepattern
            patternpos(end+1,:)=0;
        end
        arg={'patternpos',patternpos}; 
    case 'zline2'
        dim=3;
        patternpoints=length(L)*2;
        patternpos=zeros(patternpoints,3); patternpos(:,3)=[-L(1), L(1)]/2;
        if probecenter %probecenter, pattern points argument ignored if not makepattern
            patternpos(end+1,:)=0;
        end
        arg={'patternpos',patternpos}; 
    case {'octahedron'}
        dim=[1,2,3];
        patternpoints=6;
        patternpos=zeros(6,3);
        patternpos(1,1)=L/2;patternpos(2,2)=L/2; patternpos(3,1)=-L/2; patternpos(4,2)=-L/2;
        patternpos(5,3)=-L/2; patternpos(6,3)=L/2;
        if probecenter %probecenter, pattern points argument ignored if not makepattern
            patternpos(end+1,:)=0;
        end
        arg={"patternpos",patternpos}; 
    otherwise
        warning([itr.Mode.id]+ " not implemented, SimSequenceFile");                        
end
pointdwelltime=itr.patDwellTime/itr.patRepeat*1e3/patternpoints;
if probecenter
    pointdwelltime(2)=pointdwelltime(1)*patternpoints*seq.ctrDwellFactor;
end
laserpower=itr.pwrFactor;
arg2={"phasemask",phasemask, "orbitpoints",patternpoints, "orbitL",L,...
    "probecenter",probecenter,"pointdwelltime",pointdwelltime, "laserpower",...
    laserpower,"repetitions",itr.patRepeat,"pinholeorbit",pinholeorbit};
parg=horzcat(arg,arg2);

% estimators
% esth.function="est_abberior_debiased";
% esth.par={"patternpos","coefficients","patGeoFactor"};
esth.dim=dim;
switch mat2str(dim)
    case mat2str([1,2])
        % esth.dim=[1,2];
        if contains(itr.Mode.modulated,'phl')
            esth.function="est_pinholeorbit";
            esth.par={"patternpos", L, sigma_est_ph, probecenter};
        else
            esth.function="est_donutLSQ1_2D";
            esth.par={"patternpos", L, 310, 0};
        end
    case mat2str(3)
        % esth.dim=3;
        if itr.Mode.pattern=="zline" %5 points: now with 3, but make a 5 point estimaotr
            esth.function="est_zline"; % 
            esth.par={L};
        elseif itr.Mode.pattern=="zline2" %3 points
            esth.function="est_qLSQiter1D"; % 
            esth.par={L};
        end
    case mat2str([1,2,3])
        % esth.dim=[1,2,3];
        if itr.Mode.pattern=="octahedron"
            esth.par={"patternpos", L, 310, 0};
            esth.function="est_octahedron"; %
        end
end
end