%% Imaging of DNA-PAINT (blinking fluorophores + diffusive background) with Abberior scouting sequence
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path

%make abberior simulator
if ~exist('sim','var') || ~isa(sim,"SimSequencefile")
    sim=SimSequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
end
fname='Imaging_2D.json';
sim.loadsequence(fname,'PSFvectorial2D.json'); %only sequence file, then simple gauss and donut PSFs are used (fast)
% sim.makepatterns;
% sim.scoutingcoordinates=[0 0];
sim.makescoutingpattern([-80 -150; 120 100 ]) %for imaging
sim.sequence.locLimit=1000;% to avoid getting stuck with background fluorophore
%
maxtime=60*1e3; %ms, i.e. 1 min
cfrcutoff=0.5;

% make a fluorophore collection with blinking fluorophores
fc=FlCollectionBlinking;


laserpower=8;
switchpar.brightness=100*laserpower;
switchpar.toffsmlm=20*1e3; %on-switching time in ms
switchpar.photonbudget=8000;
switchpar.tonsmlm=100; % ms 
switchpar.activations=inf; %re activations
switchpar.starton=-1; %fluorophores start in random on / off state, determined by tonsmlm, toffsmlm
fc.setpar(switchpar)

%add fake NPCs
fc.add(makeNPC(pos=[0 0 0]));

% make diffusing molecules, back of the envelope calculations
% standard sequences: ~1 nM, fast sequences: 0.2 pM
% 1 M = Na/liter, 1 l= (0.1 m)^3 = (0.1 *1e6 um)^3 =1e15 um^3
% density of fluorophores (um^-3) for 1nM: 1e-9* 6e23/1e15=0.6 %density of fluorophores, about 1 / um^3
% slow sequences, bounding box: 1 um  x 1 um x 2 um: ~ 1 particle
% fast sequences: 2 um x 2 um x 2 um: ~1 particle
boundingbox=[2000 2000 2000]; %two fluorophores per 2x2x2 um
D=30; %um^2/s
fd=FlMoving;
fd.makediffusion(D,0.01,dim=3,boundarybox=boundingbox);
fd2=FlMoving;
fd2.makediffusion(D,0.01,dim=3,boundarybox=boundingbox);
fc.add({fd,fd2})

sim.fluorophores=fc;
brightnesses=[0,1]; %compare without and with background from diffusing fluorophore
titles=["DNA-PAINT: imaging strands invisible","DNA-PAINT: diffusive imaging strands"];

%
figure(300)
tiledlayout("TileSpacing","tight")

for k=1:length(brightnesses)
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
    fc.reset; %switch on all fluorophores again
    fd.brightness=switchpar.brightness*brightnesses(k);
    fd2.brightness=switchpar.brightness*brightnesses(k);
    
    out=sim.scoutingSequence(maxtime=maxtime);
    % out=sim.scoutingSequence(maxrep=repetitions(k));
    
    %plot results
    vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
    vldcfr=vld & out.loc.cfr<cfrcutoff;
    notvld=~vld & ~vldcfr;

    % subplot(2,2,k); 
    nexttile
    hold off;
    plot(sim.scoutingcoordinates(:,1),sim.scoutingcoordinates(:,2),'k*')
    hold on
    plot(out.loc.xnm(notvld),out.loc.ynm(notvld),'y.')
    plot(out.loc.xnm(vld),out.loc.ynm(vld),'m.')
    plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
    posfl=squeeze(out.fluorophores.pos(end,1:end-2,:)); %last one is diffusing
    plot(posfl(:,1),posfl(:,2),'ro')
    axis equal 
    legend('scouting','not vld', 'last itr vld','last itr vld +cfr', 'fluorophore')
    title(titles(k))
    drawnow
    plot([-90 10],[-70 -70],'k')
    text(-50,-60,"100 µm")
    ax=gca; ax.XAxis.Visible="off";ax.YAxis.Visible="off";
    xlim([-95 100])
    ylim([-85 105])
end



%% Imaging of blinking fluorophores with Abberior sequence
%make abberior simulator

photonbudget=[800, 5000];
reactivations =[0 2];
brightnesses=[50 100];
titles=["PALM","dSTORM"];

for k=1:length(photonbudget)
sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
fc=FlCollectionBlinking;
%set parameterst for caged fluorophore, PAFP or similar
% laserpower=5;
switchpar.brightness=brightnesses(k)*laserpower;
% switchpar.toffsmlm=20*1e3; %on-switching time in ms
switchpar.photonbudget=photonbudget(k);
switchpar.tonsmlm=1e8; % ms stays on, only bleached
switchpar.activations=reactivations(k); %re activations
switchpar.starton=0; %fluorophores start in random on / off state, determined by tonsmlm, toffsmlm
fc.setpar(switchpar)

%add fake NPCs
fc.addstatic(makeNPC(pos=[0 0 0]));
fc.reset;

sim.fluorophores=fc;

out=sim.scoutingSequence(maxtime=maxtime);

%plot results

vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
vldcfr=vld & out.loc.cfr<cfrcutoff;
notvld=~vld & ~vldcfr;
% subplot(2,2,k+2);
nexttile
hold off;
plot(sim.scoutingcoordinates(:,1),sim.scoutingcoordinates(:,2),'k*')
hold on
plot(out.loc.xnm(notvld),out.loc.ynm(notvld),'y.')
plot(out.loc.xnm(vld),out.loc.ynm(vld),'m.')
plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
posfl=squeeze(out.fluorophores.pos(end,:,:));
plot(posfl(:,1),posfl(:,2),'ro')
axis equal 
legend('scouting','not vld', 'last itr vld','last itr vld +cfr', 'fluorophore')
title(titles(k))
    plot([-90 10],[-70 -70],'k')
    text(-50,-60,"100 µm")
drawnow
ax=gca; ax.XAxis.Visible="off";ax.YAxis.Visible="off";
    xlim([-95 100])
    ylim([-85 105])
end



%% Fluorohore density for dSTORM and filtering
%make abberior simulator
figure(302)
tiledlayout("TileSpacing","tight")

photonbudget=5000;
reactivations =2;
brightnesses=100;
laserpower=5;

toff=[0.5, 1, 2, 5, 10, 20, 50, 100, 200]*1e3; %ms

for k=1:length(toff)
sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
fc=FlCollectionBlinking;
%set parameterst for caged fluorophore, PAFP or similar

switchpar.brightness=brightnesses*laserpower;
switchpar.toffsmlm=toff(k); %on-switching time in ms
switchpar.photonbudget=photonbudget;
switchpar.tonsmlm=1e4; % ms stays on, only bleached
switchpar.activations=reactivations; %re activations
switchpar.starton=-1; %fluorophores start in random on / off state, determined by tonsmlm, toffsmlm
fc.setpar(switchpar)

%add fake NPCs
fc.addstatic(makeNPC(pos=[0 0 0]));
fc.reset;

sim.fluorophores=fc;

out=sim.scoutingSequence(maxtime=maxtime);

%plot results

vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
vldcfr=vld & out.loc.cfr<cfrcutoff;
notvld=~vld & ~vldcfr;
% subplot(2,2,k+2);
nexttile
hold off;
plot(sim.scoutingcoordinates(:,1),sim.scoutingcoordinates(:,2),'k*')
hold on
plot(out.loc.xnm(notvld),out.loc.ynm(notvld),'y.')
plot(out.loc.xnm(vld),out.loc.ynm(vld),'m.')
plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
posfl=squeeze(out.fluorophores.pos(end,:,:));
plot(posfl(:,1),posfl(:,2),'ro')
axis equal 
legend('scouting','not vld', 'last itr vld','last itr vld +cfr', 'fluorophore')
title("toff (s): "+toff(k)/1000)
    plot([-90 10],[-70 -70],'k')
    text(-50,-60,"100 µm")
drawnow
ax=gca; ax.XAxis.Visible="off";ax.YAxis.Visible="off";
    xlim([-95 100])
    ylim([-85 105])
end