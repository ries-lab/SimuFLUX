%test simulator
%%PhaseFLUX
% clear all
% fl=MFfluorophore;
% psfphaseflux=PSFMF_PhaseFLUX;
% 
% sim=MFSimulator(fl);
%%

psfphaseflux.zerooffset=0.0;
numlocs=1000;

L=50;Lz=150;

fl.pos=[10 1 2];
fl.brightness=100;
sim.dwelltime=100;

% sim.definePattern('donut', psfdonut, makepattern='orbitscan', orbitpoints=6, probecenter=true,orbitL=L)
% sim.definePattern('vortex', psfphaseflux, psfpar="vortex", makepattern='orbitscan', orbitpoints=6, probecenter=true,orbitL=L)
sim.definePattern('x',psfphaseflux,psfpar='x',zeropos=[-1 0 1]*L/2)
sim.definePattern('y',psfphaseflux,psfpar='y',zeropos=[-1 0 1]*L/2)
sim.definePattern('z',psfphaseflux,psfpar='tophat',zeropos=[-1 0 1]*Lz/2)
sim.combinepatterns("xyz",["x","y","z"])

sim.pospattern=[0 0 0];
% sim.patternrepeat('x',1)
% sim.patternrepeat('y',5000)

xest=zeros(numlocs,3);
phot=zeros(numlocs, 3,3);
photall=zeros(numlocs,1);

for k=1:numlocs
    % photh=sim.patternrepeat(usepattern,1);
    % xest(k,:)=positionestimatedonut(photh,sim.patterns(usepattern).pos,L,psfdonut.fwhm);
    % xest(k,:)=positionestimate(photh,sim.patterns(usepattern).pos);
    photx=sim.patternrepeat("x",1);
    xest(k,:)=xest(k,:)+positionestimate1D(photx,L,1);
    photy=sim.patternrepeat("y",1);
    xest(k,:)=xest(k,:)+positionestimate1D(photy,L,2);
    photz=sim.patternrepeat("z",1);
    xest(k,:)=xest(k,:)+positionestimate1D(photz,Lz,3);
    photall(k)=sum(photx)+sum(photy)+sum(photz);
    phot(k,:,1)=photx;phot(k,:,2)=photy;phot(k,:,3)=photz;
    % phot(k)=sum(photh);
end

estimatorx=@(phot) positionestimate1D(phot,L,1);
estimatory=@(phot) positionestimate1D(phot,L,2);
estimatorz=@(phot) positionestimate1D(phot,Lz,3);

recenterh=@(x) recenter(sim,x);

seq={"x",1, estimatorx, 1, recenterh, false; "y",1, estimatory, 2, recenterh, false; "z",1, estimatorz, 3, recenterh, false};
sim.defineSequence("phaseflux3D",seq);
[xests,photons]=sim.runSequence("phaseflux3D",numlocs);
xest=xests;

% mean(phot,1)
lp=sim.calculateCRB("xyz",dim=3)/sqrt(mean(photall));

ff='%1.2f,';
disp(['mean(phot): ', num2str(mean(photall),ff),...
    ' std: ', num2str(std(xest,'omitnan'),ff),...
    ' rmse: ', num2str(rmse(xest,fl.pos,'omitnan'),ff),...
    ' pos: ', num2str(mean(xest,'omitnan'),ff),...
    ' bias: ', num2str(mean(xest-fl.pos,'omitnan'),ff),...
    ' locprec: ', num2str(sim.locprec(mean(photall),L),ff),...
    ' sqrtCRB: ', num2str(lp,ff)])

%% Donut
%%PhaseFLUX
% fl=MFfluorophore;
% psfphaseflux=PSFMF_PhaseFLUX;
% psfdonut=PSFMF_donut2D;
% sim=MFSimulator(fl);
%%
psfphaseflux.zerooffset=0.0;
numlocs=1000;

L=50;
fl.pos=[10 1 2];
fl.brightness=100;
sim.dwelltime=100;
sim.pospattern=[0 0 0];

sim.definePattern('donut', psfdonut, makepattern='orbitscan', orbitpoints=4, probecenter=true,orbitL=L)
% sim.definePattern('vortex', psfphaseflux, psfpar="vortex", makepattern='orbitscan', orbitpoints=6, probecenter=true,orbitL=L)

xest=zeros(numlocs,3);
phot=zeros(numlocs, 3,3);
photall=zeros(numlocs,1);

usepattern="donut";
for k=1:numlocs
    photh=sim.patternrepeat(usepattern,1);
    xest(k,:)=positionestimatedonut(photh,sim.patterns(usepattern).pos,L,psfdonut.fwhm);
    photall(k)=sum(photh);
    % xest(k,:)=positionestimate(photh,sim.patterns(usepattern).pos);
end

estimator=@(phot) positionestimatedonut(phot,sim.patterns(usepattern).pos,L,psfdonut.fwhm);
recenterh=@(x) recenter(sim,x);

seq={"donut",1, estimator, 1:2, recenterh, false};
sim.defineSequence("donutseq",seq);
[xests,photons]=sim.runSequence("donutseq",numlocs);
xest=xests;

% mean(phot,1)
lp=sim.calculateCRB("donut",dim=2)/sqrt(mean(photall));

ff='%1.2f,';
disp(['mean(phot): ', num2str(mean(photall),ff),...
    ' std: ', num2str(std(xest,'omitnan'),ff),...
    ' rmse: ', num2str(rmse(xest,fl.pos,'omitnan'),ff),...
    ' pos: ', num2str(mean(xest,'omitnan'),ff),...
    ' bias: ', num2str(mean(xest-fl.pos,'omitnan'),ff),...
    ' locprec: ', num2str(sim.locprec(mean(photall),L),ff),...
    ' sqrtCRB: ', num2str(lp,ff)])



function xest=positionestimatequad(photonsi,patternpos)
    pi=photonsi/sum(photonsi);
    % eq 2.63
    xest=-sum(pi'.*patternpos);
end

function xest=positionestimatedonut(photonsi,patternpos,L,fwhm)
    pi=photonsi/sum(photonsi);
    % eq 2.63
    xest=-1/(1-(L^2*log(2)/fwhm^2))*sum(pi'.*patternpos);
end

function xest=positionestimate1D(photonsi,L,coord)
xest=[0 0 0];
    % ph=photonsi/sum(photonsi);
    % eq 2.63
    xest(coord)=L/(1+sqrt(photonsi(end)/photonsi(1)))-L/2;
end

function recenter(obj,coord)
obj.pospattern=obj.pospattern+coord;
end