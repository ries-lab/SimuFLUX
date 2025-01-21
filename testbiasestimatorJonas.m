% testblinking
fl=blinkingfluorophore;
fl=staticfluorophore;
psfg=PSFMF_donut2D;

sim=MFSimulator(fl);

%%
numlocs=25;
sequencerepetitions=10;
L=151;
fl.pos=[0 0 0];
% fl.reset
fl.brightness=10000;

xp=-300:30:300;

patternrepetitions=1;
pointdwelltime=1000/patternrepetitions;
sim.posgalvo=[0 0 0];

probecenter=true;
orbitpoints=6;
sim.definePattern('vortex', psfg, psfpar="vortex", makepattern='orbitscan', orbitpoints=orbitpoints, probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime)


% betas=[1/orbitpoints 1 10 0 1000/orbitpoints];
% betas=[1 3 40 0 500*orbitpoints]
betas_4_75=[1 3 40 0 2000];
betas_6_75=[1 3 60 0 8000];
betas6_288=[1. 4  150 600 -9000];
betas=[1 3 20 0 8000];
% betas=[1 0];

clear biasx biasy
for k=1:length(xp)
    % xp(k)
    fl.pos(1)=xp(k);
    estimator=@(phot) estimators('donut2DcentermLMS',phot,sim.patterns("vortex").pos,L,310,betas);
    recenterh=@(x) recenter(sim,x);
    seq={"vortex",patternrepetitions, estimator, 1:2, recenterh, false};
    sim.defineSequence("donutseq",seq);
    out=sim.repeatSequence("donutseq",numlocs,sequencerepetitions);
    biasx(k)=mean(out.loc.xnm-out.loc.xfl);
    biasy(k)=mean(out.loc.ynm-out.loc.yfl);
    xest(k)=mean(out.loc.xnm);
    yest(k)=mean(out.loc.ynm);
end

figure(88)
hold on
plot(xp,biasx)

figure(89)
hold on
plot(xp,xest,xp,xp)

% for k=1:length(out.photch)
%     fprintf(num2str(mean(out.photch{k}),'%4.1f,'));
% end

% lp=sim.calculateCRB("vortex",dim=2)/sqrt(mean(out.loc.phot));
% 
% sim.displayresults(out, lp,L)
% 
% figure(88);hold off;plot(out.loc.time, out.loc.xnm);hold on; plot(out.loc.time,out.loc.xfl)