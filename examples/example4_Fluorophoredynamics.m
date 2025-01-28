% Fluorophpore blinking, bleaching and movement
%% flickering fluorophore and repetitions
fl=Fl_blinkbleach;
sim=Sim_Simulator(fl);
sim.fluorophores.fast_toff = 100; %off-time us
sim.fluorophores.fast_ton = 100; % on-time us

repetitions=10; %how often to repeat the pattern scan
L=75;
pointdwelltime=100/repetitions;%us
orbitorder=[1 2 3 4 5 6]; %order might matter
psf_donut=PSF_donut2D;
sim.definePattern("donut", psf_donut, makepattern="orbitscan", orbitpoints=6, ...
    probecenter=true,orbitL=L,pointdwelltime=pointdwelltime,laserpower=25,repetitions=repetitions);
sim.defineComponent("estdonut","estimator",@estimators,parameters={"donut2D",sim.patterns("donut").pos,L,360},dim=1:2);
out=sim.runSequence({"donut","estdonut"});sim.displayresults(out)

%plot std vs repetitions
allrepetitions=1:25;
stdx=zeros(length(allrepetitions),1);
stdxrel=zeros(length(allrepetitions),1);
for k=1:length(allrepetitions)
    pointdwelltime=100/allrepetitions(k);%us
    sim.definePattern("donut", psf_donut, makepattern="orbitscan", orbitpoints=6, ...
    probecenter=true,orbitL=L,pointdwelltime=pointdwelltime,laserpower=25,repetitions=allrepetitions(k));
    out=sim.runSequence({"donut","estdonut"},maxlocs=1000);
    stdx(k)=std(out.loc.xnm,'omitnan');
    crb=sim.calculateCRB("donut",dim=1:2);
    stdxrel(k)=stdx(k)/crb(1)*sqrt(mean(out.loc.phot(out.loc.phot>0)));
end


figure(134)
plot(allrepetitions,stdxrel)
xlabel('repetitions')
ylabel('std(x) /sCRB')

%% diffusing, bleaching fluorophore and tracking

