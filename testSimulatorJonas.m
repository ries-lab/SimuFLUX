% 
% fl=staticfluorophore;
% psfv=PSFMF_vectorial;
% 
% sim=MFSimulator(fl);
%%
% psfv.zerooffset=0.0;
% 
% numlocs=1000;
% 
% L=150;
% fl.pos=[50 0 0];
% fl.brightness=1000;
% laserpower=10;
% dwelltime=100;
% 
% 
% sim.definePattern('donut', psfv, phasemask="vortex", makepattern='orbitscan', orbitpoints=6, ...
%     probecenter=0,orbitL=L,pointdwelltime=dwelltime,laserpower=laserpower,repetitions=2)
% sim.defineComponent("estdonut","estimator",@estimators,parameters={"donut2D",sim.patterns(usepattern).pos,L,fwhm},dim=1:2);
% sim.defineComponent("recenter","centerer",@simplerecenter,dim=1:2);
% 
% fwhm=360;
% 
% seq={"donut","estdonut"};
% sim.defineSequence("donutseq",seq);
% out=sim.runSequence("donutseq",numlocs);
% 
% 
% % mean(phot,1)
% lp=sim.calculateCRB("donut",dim=2)/sqrt(mean(out.loc.phot));
% sim.displayresults(out,lp,L)
% 

% adsfdf
% ff='%1.2f,';
% disp(['mean(phot): ', num2str(mean(photall),ff),...
%     ' std: ', num2str(std(xest,'omitnan'),ff),...
%     ' rmse: ', num2str(rmse(xest,fl.pos,'omitnan'),ff),...
%     ' pos: ', num2str(mean(xest,'omitnan'),ff),...
%     ' bias: ', num2str(mean(xest-fl.pos,'omitnan'),ff),...
%     ' locprec: ', num2str(sim.locprec(mean(photall),L),ff),...
%     ' sqrtCRB: ', num2str(lp,ff)])


% 
% function xest=positionestimatequad(photonsi,patternpos)
%     pi=photonsi/sum(photonsi);
%     % eq 2.63
%     xest=-sum(pi'.*patternpos);
% end
% 
% function xest=positionestimatedonut(photonsi,patternpos,L,fwhm)
%     pi=photonsi/sum(photonsi);
%     % eq 2.63
%     xest=-1/(1-(L^2*log(2)/fwhm^2))*sum(pi'.*patternpos);
% end
% 
% function xest=positionestimate1D(photonsi,L,coord)
% xest=[0 0 0];
%     % ph=photonsi/sum(photonsi);
%     % eq 2.63
%     xest(coord)=L/(1+sqrt(photonsi(end)/photonsi(1)))-L/2;
% end

% function recenter(obj,coord)
% obj.pospattern=obj.pospattern+coord;
% end

%%

psfv.zerooffset=0.0;
numlocs=1000;
dwelltime=100;

L=100;Lz=300;
laserpower=10;
sim.fluorophores.pos=[30 3 3];
fwhm=360;
% fl.brightness=100;
% sim.dwelltime=100;

% sim.definePattern('donut', psfv, phasemask="vortex", makepattern='orbitscan', orbitpoints=6, ...
%     probecenter=0,orbitL=L,pointdwelltime=dwelltime,laserpower=laserpower,repetitions=2)
% sim.definePattern('donut', psfdonut, makepattern='orbitscan', orbitpoints=6, probecenter=true,orbitL=L)
% sim.definePattern('vortex', psfphaseflux, phasemask="vortex", makepattern='orbitscan', orbitpoints=6, probecenter=true,orbitL=L)
sim.definePattern('x',psfv,phasemask='halfmoonx',zeropos=[-1 0 1]*L/2,pointdwelltime=dwelltime,laserpower=laserpower)
sim.definePattern('y',psfv,phasemask='halfmoony',zeropos=[-1 0 1]*L/2,pointdwelltime=dwelltime,laserpower=laserpower)
sim.definePattern('z',psfv,phasemask='tophat',zeropos=[-1 0 1]*Lz/2,pointdwelltime=dwelltime,laserpower=laserpower)


sim.defineComponent("estx","estimator",@estimators,parameters={"donut1D",[-1; 0 ;1]*L/2,fwhm},dim=1);
sim.defineComponent("esty","estimator",@estimators,parameters={"donut1D",[-1; 0 ;1]*L/2,fwhm},dim=2);
sim.defineComponent("estz","estimator",@estimators,parameters={"donut1D",[-1; 0 ;1]*Lz/2,fwhm},dim=3);
sim.defineComponent("posupdate","positionupdater",@simple_recenter,dim=1:2);

seq={'x','estx','y','esty','z','estz',"posupdate"};

% recenterh=@(x) recenter(sim,x);

% seq={"x",1, estimatorx, 1, recenterh, false; "y",1, estimatory, 2, recenterh, false; "z",1, estimatorz, 3, recenterh, false};
% sim.defineSequence("phaseflux3D",seq);
out=sim.runSequence(seq,maxlocs=numlocs,repetitions=2);


% lp=sim.calculateCRB("donut",dim=2)/sqrt(mean(out.loc.phot));
sim.displayresults(out,[],L)
% %% Donut
% %%PhaseFLUX
% fl=MFfluorophore;
% psfv=PSFMF_vectorial;
% % psfdonut=PSFMF_donut2D;
% sim=MFSimulator(fl);
% psfv.setpar('beadradius',00*1e-9); %sys: in m
