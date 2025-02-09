%% Imaging of DNA-PAINT (static fluorophores + diffusive background) 

if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; 
end
psf_vec.setpinhole("AU",1);
L=75;
D=30; %um^2/s
numberOfLocalizations=1000;

fd=FlMoving;
fd.makediffusion(D,0.1,dim=3,boundarybox=[500 500 1000]);
fd.brightness=10000;
fs=FlStatic;
fs.pos=[10 0 0];
fs.brightness=10000;
fc=FlCollectionBlinking;
fc.add({fs, fd});
sim1=Simulator(fc);

sim1.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=4, ...
    probecenter=false,orbitL=L,pointdwelltime=0.1,laserpower=1,repetitions=1)
fwhm=360;% size of the donut, needed for proper estimation. 
sim1.defineComponent("estdonut","estimator",@est_donut2d,parameters={sim1.patterns("donut").pos,L,fwhm},dim=1:2);

seq={"donut","estdonut"};
out=sim1.runSequence(seq,"maxlocs",numberOfLocalizations);
sim1.summarize_results(out); %display summary of simulation

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
laserpower=2;
switchpar.brightness=100*laserpower;
switchpar.toffsmlm=2*1e3; %on-switching time in ms
switchpar.photonbudget=5000;
switchpar.tonsmlm=100; % ms 
switchpar.activations=inf; %re activations
switchpar.starton=0; %fluorophores start in random on / off state, determined by tonsmlm, toffsmlm
fc.setpar(switchpar)

%add fake NPCs
fc.add(posnpc);

% make diffusing molecules, back of the envelope calculations
% standard sequences: ~1 nM, fast sequences: 0.2 pM
% 1 M = Na/liter, 1 l= (0.1 m)^3 = (0.1 *1e6 um)^3 =1e15 um^3
% c1nM= 1e-9* 6e23/1e15 %density of fluorophores, about 1 / um^3
% slow sequences, bounding box: 1 um  x 1 um x 2 um: ~ 1 particle
% fast sequences: 2 um x 2 um x um: ~1 particle
D=30; %um^2/s
fd=FlMoving;

fd.makediffusion(D,0.01,dim=3,boundarybox=[500 500 1000]);
fc.add(fd)

sim.fluorophores=fc;
maxtime=10*1e3; %10 seconds
% repetitions=[500,100];
brightnesses=[0,1]; %compare without and with background from diffusing fluorophore
% sim.bgcSenseValue=0;
for k=1:length(brightnesses)
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];sim.time=0;
    fd.brightness=switchpar.brightness*brightnesses(k);
    
    out=sim.scoutingSequence(maxtime=maxtime);
    % out=sim.scoutingSequence(maxrep=repetitions(k));
    
    %plot results
    vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
    vldcfr=vld & out.loc.cfr<0.1;
    figure(89+k); hold off;
    plot(sim.scoutingcoordinates(:,1),sim.scoutingcoordinates(:,2),'k+')
    hold on
    plot(out.loc.xnm(vld),out.loc.ynm(vld),'r.')
    plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
    posfl=squeeze(out.fluorophores.pos(end,1:end-1,:)); %last one is diffusing
    plot(posfl(:,1),posfl(:,2),'mo')
    axis equal
    legend('scouting','last itr vld','last itr vld +cfr', 'fluorophore')
end