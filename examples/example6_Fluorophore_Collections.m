% Simulate several fluorophores in a fluorophore collection
%% Imaging of blinking fluorophores with Abberior sequence
%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"Sim_sequencefile")
    sim=Sim_sequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
end
fname='Imaging_2D.json';
sim.loadsequence(fname); %only sequence file, then simple gauss and donut PSFs are used (fast)
sim.makepatterns;
sim.makescoutingpattern([-100 -100; 400 250 ]) %for imaging

% make fake NPCs
clear posnpc;
R=50;dphi=pi/4;
phi=0:dphi:2*pi-dphi; posnpc(:,1)=R*cos(phi); posnpc(:,2)= R*sin(phi);

% make a fluorophore collection with blinking fluorophores
fc=Flcollection_blinking;

%set parameterst for caged fluorophore, PAFP or similar
laserpower=5;
switchpar.brightness=100*laserpower;
switchpar.toffsmlm=3*1e6; %on-switching time in ns
switchpar.photonbudget=5000;
switchpar.tonsmlm=1e8; % ns stays on, only bleached
switchpar.activations=5; %re activations
switchpar.starton=0; %fluorophores start in random on / off state, determined by tonsmlm, toffsmlm
fc.setpar(switchpar)

%add fake NPCs
fc.add(posnpc);
% add more NPCs at positions dpos
dpos=[250 50 ]; fc.addstatic(posnpc+dpos);
dpos=[50 150 ]; fc.addstatic(posnpc+dpos);

sim.fluorophores=fc;
out=sim.scoutingSequence(maxrep=200);

%plot results
vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
vldcfr=vld & out.loc.cfr<0.1;
figure(89); hold off;
plot(sim.scoutingcoordinates(:,1),sim.scoutingcoordinates(:,2),'k+')
hold on
plot(out.loc.xnm(vld),out.loc.ynm(vld),'r.')
plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
posfl=vertcat(fc.flall.pos);
plot(posfl(:,1),posfl(:,2),'mo')
axis equal
legend('scouting','last itr vld','last itr vld +cfr', 'fluorophore')