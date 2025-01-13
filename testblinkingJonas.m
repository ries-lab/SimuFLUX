% testblinking
% fl=blinkingfluorophore;
fl=staticfluorophore;
psfdonut=PSFMF_donut2D;
sim=MFSimulator(fl);
numlocs=12;
%%
L=50;
fl.pos=[0 0 0];
% fl.toff=100;
fl.brightness=100000;

% fl.starton=false;
% fl.photonbudget=inf; 


repetitions=1;
sim.dwelltime=10/repetitions;
sim.pospattern=[0 0 0];

fl.pos={@(t) -0.05*t+2,@(t) 0};
fl.posmode='function';

% fl.makediffusion(.02,sim.dwelltime,dim=2)


sim.definePattern('donut', psfdonut, makepattern='orbitscan', orbitpoints=4, probecenter=false,orbitL=L,orbitorder=[1 2 4 3])
estimator=@(phot) estimators('donut2D',phot,sim.patterns("donut").pos,L,psfdonut.fwhm);
recenterh=@(x) recenter(sim,x);
seq={"donut",repetitions, estimator, 1:2, recenterh, true};
sim.defineSequence("donutseq",seq);
out=sim.runSequence("donutseq",numlocs);

for k=1:length(photons)
    fprintf(num2str(mean(photons{k}),'%4.1f,'));
end
% lp=sim.calculateCRB("donut",dim=2)/sqrt(mean(photall));

sim.displayresults(out, lp,L)

figure(88);hold off;plot(out.time, out.xest(:,1));hold on; plot(out.time,out.flpos(:,1))