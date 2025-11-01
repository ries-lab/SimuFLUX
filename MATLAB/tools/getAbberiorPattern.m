function [parg,esth]=getAbberiorPattern(itr,seq)
if contains(itr.Mode.modulated,'phl')&& contains(itr.Mode.pattern,'hexagon') %pinhole orbit: for now Gauss, later implement as pattern, remove here
    phasemask="flat";
elseif contains(itr.Mode.epsf,'focFldRing')
    phasemask="tophat";
elseif contains(itr.Mode.epsf,'focFldVortex')
    phasemask="vortex";
end

% scan patterns
if itr.ccrLimit==-1
    probecenter=false;
else
    probecenter=true;
end
L=itr.patGeoFactor*360; %nm
arg2D={'makepattern','orbitscan'};
dim="xy";
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
        dim="z";
        patternpoints=length(L)*2;
        patternpos=zeros(patternpoints,3); patternpos(:,3)=[-L(1), L(1), -L(2), L(2)]/2;
        if probecenter %probecenter, pattern points argument ignored if not makepattern
            patternpos(end+1,:)=0;
        end
        arg={'patternpos',patternpos}; 
    case 'zline2'
        dim="z";
        patternpoints=length(L)*2;
        patternpos=zeros(patternpoints,3); patternpos(:,3)=[-L(1), L(1)]/2;
        if probecenter %probecenter, pattern points argument ignored if not makepattern
            patternpos(end+1,:)=0;
        end
        arg={'patternpos',patternpos}; 
    case {'octahedron'}
        dim="xyz";
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
    "probecenter",probecenter,"pointdwelltime",pointdwelltime, "laserpower",laserpower,"repetitions",itr.patRepeat };
parg=horzcat(arg,arg2);

% estimators
switch dim
    case"xy" 
        esth.dim=[1,2];
        if contains(itr.Mode.modulated,'phl')
            esth.function="est_GaussLSQ1_2D";
            esth.par={"patternpos", L, 120, probecenter};
        else
            esth.function="est_donutLSQ1_2D";
            esth.par={"patternpos", L, 310, 0};
        end
    case "z"
        esth.dim=3;
        if itr.Mode.pattern=="zline" %5 points: now with 3, but make a 5 point estimaotr
            esth.function="est_zline"; % 
            esth.par={L};
        elseif itr.Mode.pattern=="zline2" %3 points
            esth.function="est_qLSQiter1D"; % 
            esth.par={L};
        end
    case "xyz"
        esth.dim=[1,2,3];
        if itr.Mode.pattern=="octahedron"
            esth.par={"patternpos", L, 310, 0};
            esth.function="est_octahedron"; %
        end
end
end