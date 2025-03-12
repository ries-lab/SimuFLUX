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
sim1=Simulator(fluorophores=fc);

sim1.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=4, ...
    probecenter=true,orbitL=L,pointdwelltime=0.1,laserpower=1,repetitions=1)
fwhm=360;% size of the donut, needed for proper estimation. 
sim1.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,true},dim=1:2);

seq={"donut","estdonut"};
out=sim1.runSequence(seq,"maxlocs",numberOfLocalizations);
sim1.summarize_results(out); %display summary of simulation

%% Imaging of DNA-PAINT (blinking fluorophores + diffusive background) with Abberior scouting sequence
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path

%make abberior simulator
if ~exist('sim2','var') || ~isa(sim2,"SimSequencefile")
    sim2=SimSequencefile;
else
    sim2.posgalvo=[0 0 0];sim2.posEOD=[0 0 0];sim2.time=0;
end
fname='Imaging_2D.json';
sim2.loadsequence(fname,'PSFvectorial2D.json'); %only sequence file, then simple gauss and donut PSFs are used (fast)
% sim.makepatterns;
% sim.scoutingcoordinates=[0 0];
sim2.makescoutingpattern([-100 -150; 100 100 ]) %for imaging
sim2.sequence.locLimit=100;% to avoid getting stuck with background fluorophore
%

% make a fluorophore collection with blinking fluorophores
fc=FlCollectionBlinking;

%set parameterst for caged fluorophore, PAFP or similar
laserpower=8;
switchpar.brightness=100*laserpower;
switchpar.toffsmlm=15*1e3; %on-switching time in ms
switchpar.photonbudget=5000;
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
D=30; %um^2/s
fd=FlMoving;

fd.makediffusion(D,0.01,dim=3,boundarybox=[500 500 1000]);
fc.add(fd)

sim2.fluorophores=fc;
maxtime=30*1e3; %10 seconds
cfrcutoff=0.5;
brightnesses=[0,1]; %compare without and with background from diffusing fluorophore
titles=["imaging strands invisible","diffusive imaging strands"];

for k=1:length(brightnesses)
    sim2.posgalvo=[0 0 0];sim2.posEOD=[0 0 0];sim2.time=0;
    fc.reset; %switch on all fluorophores again
    fd.brightness=switchpar.brightness*brightnesses(k);
    
    out=sim2.scoutingSequence(maxtime=maxtime);
    % out=sim.scoutingSequence(maxrep=repetitions(k));
    
    %plot results
    vld=out.loc.vld==1 & out.loc.itr==max(out.loc.itr) ;
    vldcfr=vld & out.loc.cfr<cfrcutoff;
    figure(290+k); hold off;
    plot(sim2.scoutingcoordinates(:,1),sim2.scoutingcoordinates(:,2),'k*')
    hold on
    plot(out.loc.xnm(vld),out.loc.ynm(vld),'m.')
    plot(out.loc.xnm(vldcfr),out.loc.ynm(vldcfr),'bx')
    posfl=squeeze(out.fluorophores.pos(end,1:end-1,:)); %last one is diffusing
    plot(posfl(:,1),posfl(:,2),'ro')
    axis equal
    legend('scouting','last itr vld','last itr vld +cfr', 'fluorophore')
    title(titles(k))
    drawnow
end