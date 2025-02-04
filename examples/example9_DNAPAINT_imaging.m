%% Imaging of DNA-PAINT (blinking fluorophores + diffusive background) with Abberior sequence
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path

%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"SimSequencefile")
    sim=SimSequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
end
fname='Imaging_2D.json';
sim.loadsequence(fname,'PSFvectorial2D.json'); %only sequence file, then simple gauss and donut PSFs are used (fast)
sim.makepatterns;
sim.scoutingcoordinates=[0 0];
sim.makescoutingpattern([-1 -1; 1 1 ]*200) %for imaging
sim.sequence.locLimit=100;% to avoid getting stuck with background fluorophore
%
% make fake NPCs
clear posnpc;
R=50;dphi=pi/4;
phi=0:dphi:2*pi-dphi; posnpc(:,1)=R*cos(phi); posnpc(:,2)= R*sin(phi);

% make a fluorophore collection with blinking fluorophores
fc=FlCollectionBlinking;

%set parameterst for caged fluorophore, PAFP or similar
laserpower=5;
switchpar.brightness=100*laserpower;
switchpar.toffsmlm=5*1e3; %on-switching time in ms
switchpar.photonbudget=5000;
switchpar.tonsmlm=100; % ms stays on, only bleached
switchpar.activations=5; %re activations
switchpar.starton=0; %fluorophores start in random on / off state, determined by tonsmlm, toffsmlm
fc.setpar(switchpar)

%add fake NPCs
fc.add(posnpc);

% make diffusing molecules
% standard sequences: ~1 nM
% fast sequences: 0.2 pM
% 1 M = Na/liter 
% 1 l= (0.1 m)^3 = (0.1 *1e6 um)^3 =1e15 um^3
% c1nM= 1e-9* 6e23/1e15
% bounding box: 1 um  x 1 um x 2 um: ~ 1 particle
% fast sequences 2 um x 2 um x um: ~1 particle
D=10; %um^2/s
fd=FlMoving;
fd.brightness=switchpar.brightness;
fd.makediffusion(D,0.01,dim=3,boundarybox=[.5 0.5 1]);
fc.add(fd)

sim.fluorophores=fc;

out=sim.scoutingSequence(maxrep=10);

%plot results
vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
vldcfr=vld & out.loc.cfr<0.1;
figure(89); hold off;
plot(sim.scoutingcoordinates(:,1),sim.scoutingcoordinates(:,2),'k+')
hold on
plot(out.loc.xnm(vld),out.loc.ynm(vld),'r.')
plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
posfl=out.fluorophores.pos(end,:,2:end);
plot(posfl(:,1),posfl(:,2),'mo')
posdiff=squeeze(out.fluorophores.pos(:,end,:));
plot(posdiff(:,1),posdiff(:,2))
axis equal
legend('scouting','last itr vld','last itr vld +cfr', 'fluorophore')