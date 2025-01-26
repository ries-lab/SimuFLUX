% % % testblinking
% fl=blinkingfluorophore;
% % fl=staticfluorophore;
% psfg=PSFMF_vectorial;
% 
% sim=MFSimulator(fl);

%%
numlocs=12;
sequencerepetitions=100;
L=75;
fl.pos=[10 0 0];
fl.reset
% fl.toff=0;
fl.brightness=10000;

psfg.setpinhole("AU",1,offset=[50 0]);

% fl.starton=false;
% fl.photonbudget=inf; 


patternrepetitions=1;
pointdwelltime=1000/patternrepetitions;
sim.posgalvo=[0 0 0];
% 
% fl.pos={@(t) -0.05*t+2,@(t) 0};
% fl.posmode='function';

% fl.makesteps(16,10,4,angle=45,numpoints=100);

% fl.makediffusion(.02,sim.dwelltime,dim=2)
probecenter=false;
psfg.parameters.sys.maskshift=[0.1,0.1];
sim.definePattern('vortex', psfg, phasemask="vortex", makepattern='orbitscan', orbitpoints=6, probecenter=true,orbitL=L,pointdwelltime=pointdwelltime)

% sim.definePattern('donut', psfdonut, makepattern='orbitscan', orbitpoints=6, probecenter=probecenter,...
    % orbitL=L,pointdwelltime=pointdwelltime)
estimator=@(phot) estimators('donut2D',phot,sim.patterns("vortex").pos,L,360,probecenter);
recenterh=@(x) recenter(sim,x);
seq={"vortex",patternrepetitions, estimator, 1:2, recenterh, false};
sim.defineSequence("donutseq",seq);
out=sim.repeatSequence("donutseq",numlocs,sequencerepetitions);

% for k=1:length(out.photch)
%     fprintf(num2str(mean(out.photch{k}),'%4.1f,'));
% end

lp=sim.calculateCRB("vortex",dim=2)/sqrt(mean(out.loc.phot));

sim.displayresults(out, lp,L)

figure(88);hold off;plot(out.loc.time, out.loc.xnm);hold on; plot(out.loc.time,out.loc.xfl)