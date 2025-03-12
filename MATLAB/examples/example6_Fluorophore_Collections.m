% Simulate several fluorophores in a fluorophore collection
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
%% two fluorophores
fl1=FlStatic;
fl2=FlStatic;
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);

fc=FlCollection;
fc.add({fl1, fl2})

sim=Simulator(fluorophores=fc); %make a simulator and attach fluorophore

numberOfLocalizations=1000;
orbitpoints=4;
probecenter=true;
L=75;
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,laserpower=100)
% sim.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim.patterns("donut").pos,L,360},dim=1:2);
sim.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);


seq={"donut","estdonut"};
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
stats=sim.summarize_results(out);

x=0:20:750;

figure(260)
tiledlayout(1,2,"TileSpacing","tight")
nexttile
sim.scan_fov(seq,x,dimplot=1,dimscan=1,fluorophorenumber=2,ax1=["std","rmse","sCRB","bias"],title="scan second fluorophore in x", tag="x");



% now z-dependence
fl1.pos=[0 0 0]; fl2.pos=[100 0 0];
z=0:50:1000;
nexttile
sim.scan_fov(seq,z,dimplot=1,dimscan=3,fluorophorenumber=2,ax1=["std","rmse","sCRB","bias"],title="scan second fluorophore in z",tag="z");



%% Imaging of blinking fluorophores with Abberior sequence
%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"SimSequencefile")
    sim=SimSequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;sim.background=0;
end
fname='Imaging_2D.json';
sim.loadsequence(fname); %only sequence file, then simple gauss and donut PSFs are used (fast)
% sim.makepatterns;
sim.makescoutingpattern([-100 -100; 400 250 ]) %for imaging

% make fake NPCs

% make a fluorophore collection with blinking fluorophores


photonbudget=[800, 5000];
reactivations =[0 3];
titles=["PALM","dSTORM"];

for k=1:length(photonbudget)
fc=FlCollectionBlinking;
%set parameterst for caged fluorophore, PAFP or similar
laserpower=5;
switchpar.brightness=100*laserpower;
switchpar.toffsmlm=10*1e3; %on-switching time in ms
switchpar.photonbudget=photonbudget(k);
switchpar.tonsmlm=1e8; % ms stays on, only bleached
switchpar.activations=reactivations(k); %re activations
switchpar.starton=0; %fluorophores start in random on / off state, determined by tonsmlm, toffsmlm
fc.setpar(switchpar)

%add fake NPCs
fc.addstatic(makeNPC(pos=[0 0 0]));
% add more NPCs at positions dpos
fc.addstatic(makeNPC(pos=[250 50 0]));
fc.addstatic(makeNPC(pos=[50 150 0]));

sim.fluorophores=fc;
out=sim.scoutingSequence(maxrep=5000);

%plot results

vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
vldcfr=vld & out.loc.cfr<0.1;
notvld=~vld & ~vldcfr;
figure(262+k); hold off;
plot(sim.scoutingcoordinates(:,1),sim.scoutingcoordinates(:,2),'k+')
hold on
plot(out.loc.xnm(notvld),out.loc.ynm(notvld),'c.')
plot(out.loc.xnm(vld),out.loc.ynm(vld),'m.')
plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
posfl=squeeze(out.fluorophores.pos(end,:,:));
plot(posfl(:,1),posfl(:,2),'ro')
axis equal
legend('scouting','not vld', 'last itr vld','last itr vld +cfr', 'fluorophore')
title(titles(k))
drawnow
end