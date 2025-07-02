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
% orbitorder=[1 2 3 4 5 6 7]; %order might matter. If you change this, you cannot use est_qLSQiter2D.
orbitpoints=4;
probecenter=true;
psf_donut=PsfDonut2D;
sim.definePattern("donut", psf_donut, makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=30,...
    repetitions=repetitions);
% sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
sim.defineComponent("estsq","estimator",@est_qLSQiter2D,parameters={L,probecenter},dim=1:2);
out=sim.runSequence({"donut","estsq"});sim.summarize_results(out);

%% plot std vs repetitions
%Fig. 2
allrepetitions=1:1:25;
laserpowers=[10 100 100000]/0.88;


clear results
for l=1:length(laserpowers)
    for k=1:length(allrepetitions)
        pointdwelltime=pointdwelltimerep/allrepetitions(k);%us
        sim.definePattern("donut4", psf_donut, makepattern="orbitscan", orbitpoints=orbitpoints, ...
        probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpowers(l),repetitions=allrepetitions(k));
        out=sim.runSequence({"donut4","estsq"},maxlocs=3000);
        bright=out.loc.phot>quantile(out.loc.phot,0.1); %filter out localizations that are too dim, outliers from fluorophores that are mostly off
        converged=abs(out.loc.xnm)<L & abs(out.loc.ynm)<L; %estimator is artificially bound by +/-L. Here we remove localizations where the estimator converged to this limit.
        sr=sim.summarize_results(out,filter=bright&converged,display=false);
        results(l).stdx(k)=sr.std(1);
        results(l).stdy(k)=sr.std(2);
        % crb=sr.sCRB(1);
        results(l).stdxrel(k)=stdx(k)/sr.sCRB(1);
        results(l).stdyrel(k)=stdy(k)/sr.sCRB(2);
        results(l).biasx(k)=sr.bias(1);
        results(l).biasy(k)=sr.bias(2);
        results(l).crb(k)=sr.sCRB(1);
        results(l).phot(k)=sr.phot;
    end
end
%%
figure(240)
hold off
plot(allrepetitions,results(1).stdx,'rx',allrepetitions,results(1).crb,'r--')
hold on
plot(allrepetitions,results(2).stdx,'bo',allrepetitions,results(2).crb,'b--')
plot(allrepetitions,results(3).stdx,'k',allrepetitions, sqrt(results(3).stdx.^2+results(1).crb.^2),'r',allrepetitions, sqrt(results(3).stdx.^2+results(2).crb.^2),'b')

xlabel('repetitions')
ylabel('std (nm)')
title("flickering fluorophores, mitigation by repetitive scanning")
legend("N="+round(results(1).phot(1)), "CRB", "N="+round(results(2).phot(1)), "CRB", "N high", "prediction", "prediction")

