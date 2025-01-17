% testblinking
fl=blinkingfluorophore;
% fl=staticfluorophore;
psfdonut=PSFMF_donut2D;
sim=MFSimulator(fl);
numlocs=12;
%%
L=150;
fl.pos=[0 0 0];
fl.toff=0;
fl.brightness=1000;

% fl.starton=false;
% fl.photonbudget=inf; 


repetitions=1;
pointdwelltime=100/repetitions;
sim.pospattern=[0 0 0];
% 
% fl.pos={@(t) -0.05*t+2,@(t) 0};
% fl.posmode='function';

% fl.makesteps(16,10,4,angle=45,numpoints=100);

% fl.makediffusion(.02,sim.dwelltime,dim=2)


sim.definePattern('donut', psfdonut, makepattern='orbitscan', orbitpoints=4, probecenter=false,...
    orbitL=L,orbitorder=[1 2 4 3],pointdwelltime=pointdwelltime)
estimator=@(phot) estimators('donut2D',phot,sim.patterns("donut").pos,L,psfdonut.fwhm);
recenterh=@(x) recenter(sim,x);
seq={"donut",repetitions, estimator, 1:2, recenterh, true};
sim.defineSequence("donutseq",seq);
out=sim.runSequence("donutseq",numlocs);

% for k=1:length(out.photch)
%     fprintf(num2str(mean(out.photch{k}),'%4.1f,'));
% end

lp=sim.calculateCRB("donut",dim=2)/sqrt(mean(out.loc.phot));

sim.displayresults(out, lp,L)

figure(88);hold off;plot(out.loc.time, out.loc.xnm);hold on; plot(out.loc.time,out.loc.xfl)