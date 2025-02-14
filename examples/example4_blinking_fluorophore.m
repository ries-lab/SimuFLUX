% Fluorophpore blinking, bleaching and movement
%% flickering fluorophore and repetitions
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
fl=FlBlinkBleach;
sim=Simulator(fluorophores=fl);
sim.fluorophores.fast_toff = .1; %off-time ms
sim.fluorophores.fast_ton = .1; % on-time ms

repetitions=1; %how often to repeat the pattern scan
L=75;
pointdwelltimerep=0.1; %ms, for all repetitions
pointdwelltime=pointdwelltimerep/repetitions;%ms
orbitorder=[1 2 3 4 5 6 7]; %order might matter. If you change this, you cannot use est_quad2Diter.
orbitpoints=6;
probecenter=true;
psf_donut=PsfDonut2D;
sim.definePattern("donut", psf_donut, makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=100,...
    repetitions=repetitions,orbitorder=orbitorder);
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
out=sim.runSequence({"donut","estdonut"});sim.summarize_results(out);

%% plot std vs repetitions
allrepetitions=1:1:25;
stdx=zeros(length(allrepetitions),1);stdy=stdx;stdxrel=stdx;stdyrel=stdy;biasx=stdy;biasy=stdy;
sim.defineComponent("estsq","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
for k=1:length(allrepetitions)
    pointdwelltime=pointdwelltimerep/allrepetitions(k);%us
    sim.definePattern("donut4", psf_donut, makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=100,repetitions=allrepetitions(k));
    out=sim.runSequence({"donut4","estsq"},maxlocs=3000);
    bright=out.loc.phot>quantile(out.loc.phot,0.05); %filter out localizations that are too dim, outliers from fluorophores that are mostly off
    sr=sim.summarize_results(out,filter=bright,display=false);
    stdx(k)=sr.std(1);
    stdy(k)=sr.std(2);
    % crb=sr.sCRB(1);
    stdxrel(k)=stdx(k)/sr.sCRB(1);
    stdyrel(k)=stdy(k)/sr.sCRB(1);
    biasx(k)=sr.bias(1);
    biasy(k)=sr.bias(2);
end
figure(240)
plot(allrepetitions,stdxrel,allrepetitions,stdyrel)
xlabel('repetitions')
ylabel('std/sCRB')

