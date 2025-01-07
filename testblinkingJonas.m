% testblinking
fl=blinkingfluorophore;
fl=staticfluorophore;
psfdonut=PSFMF_donut2D;
sim=MFSimulator(fl);
numlocs=5000;
%%
L=50;
fl.pos=[0 0 0];
% fl.toff=100;
fl.brightness=1000;
% fl.toff=0;
% fl.starton=true;
% fl.photonbudget=inf; 
repetitions=1;
sim.dwelltime=10/repetitions;
sim.pospattern=[0 0 0];

sim.definePattern('donut', psfdonut, makepattern='orbitscan', orbitpoints=4, probecenter=true,orbitL=L)
estimator=@(phot) estimators('donut2D',phot,sim.patterns("donut").pos,L,psfdonut.fwhm);
recenterh=@(x) recenter(sim,x);
seq={"donut",repetitions, estimator, 1:2, recenterh, false};
sim.defineSequence("donutseq",seq);
[xest,photons,photall]=sim.runSequence("donutseq",numlocs);

for k=1:length(photons)
    fprintf(num2str(mean(photons{k}),'%4.0f,'));
end
lp=sim.calculateCRB("donut",dim=2)/sqrt(mean(photall));
sim.displayresults(xest, photall, lp,L)

