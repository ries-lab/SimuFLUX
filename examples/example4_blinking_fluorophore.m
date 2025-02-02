% Fluorophpore blinking, bleaching and movement
%% flickering fluorophore and repetitions
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
fl=FlBlinkBleach;
sim=Simulator(fl);
sim.fluorophores.fast_toff = .1; %off-time ms
sim.fluorophores.fast_ton = .1; % on-time ms

repetitions=10; %how often to repeat the pattern scan
L=75;
pointdwelltimerep=0.1; %ms, for all repetitions
pointdwelltime=pointdwelltimerep/repetitions;%ms
orbitorder=[1 2 3 4 5 6]; %order might matter
orbitpoints=6;
probecenter=true;
psf_donut=PsfDonut2D;
sim.definePattern("donut", psf_donut, makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=25,repetitions=repetitions);
sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
out=sim.runSequence({"donut","estdonut"});sim.summarize_results(out);

%plot std vs repetitions
allrepetitions=1:5:25;
stdx=zeros(length(allrepetitions),1);stdy=stdx;stdxrel=stdx;stdyrel=stdy;biasx=stdy;biasy=stdy;
for k=1:length(allrepetitions)
    pointdwelltime=pointdwelltimerep/allrepetitions(k);%us
    sim.definePattern("donut", psf_donut, makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=25,repetitions=allrepetitions(k));
    out=sim.runSequence({"donut","estdonut"},maxlocs=1000);
    sr=sim.summarize_results(out);
    stdx(k)=sr.std(1);
    stdy(k)=sr.std(2);
    % crb=sr.sCRB(1);
    stdxrel(k)=stdx(k)/sr.sCRB(1);
    stdyrel(k)=stdy(k)/sr.sCRB(1);
    biasx(k)=sr.bias(1);
    biasy(k)=sr.bias(2);

    % stdx(k)=std(out.loc.xnm,'omitnan');
    % stdy(k)=std(out.loc.ynm,'omitnan');
    % crb=sim.calculateCRB("donut",dim=1:2);
    % stdxrel(k)=stdx(k)/crb(1)*sqrt(mean(out.loc.phot(out.loc.phot>0)));
    % stdyrel(k)=stdy(k)/crb(2)*sqrt(mean(out.loc.phot(out.loc.phot>0)));
    % biasx(k)=mean(out.loc.xnm-out.loc.xfl1,'omitnan');
    % biasy(k)=mean(out.loc.ynm-out.loc.yfl1,'omitnan');
end
figure(134)
plot(allrepetitions,stdxrel,allrepetitions,stdyrel)
xlabel('repetitions')
ylabel('std/sCRB')

