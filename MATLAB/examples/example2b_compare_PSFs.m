%% Vectorial PSF
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.setpinhole("AU",1);
psf_vec.zerooffset=0.000; %true zero
sys_aberr.Zr(1,1)=4;sys_aberr.Zr(1,2)=0;sys_aberr.Zr(1,3)=0.0; %spherical aberrations 
psf_vec.setpar(sys_aberr)

fl=FlStatic;
fl.pos=[10 0 0];
sim=Simulator(fluorophores=fl); %make a simulator and attach fluorophore

maxlocs=5000;

%define scan pattern
L=75; %size of scan pattern
zeroposx=[-1;1;0]*L/2;
probecenter=true; %should we also probe the center?
orbitpoints=6;
laserpowerdonut=10; %relative, increases brightness
laserpowertophat=80;
laserpowerpf=8;
pointdwelltime=0.1; % ms, measurement time in each point
repetitions=1; %how often to repeat the pattern scan
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpowerdonut,repetitions=repetitions)
sim.definePattern("tophat_xy", psf_vec, phasemask="tophat", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpowertophat,repetitions=repetitions,dim=1:2);
sim.definePattern("pf_x", psf_vec, phasemask="halfmoonx", zeropos=zeroposx,...
    pointdwelltime=pointdwelltime,laserpower=laserpowerpf,repetitions=repetitions,dim=1);
sim.definePattern("pf_y", psf_vec, phasemask="halfmoony", zeropos=zeroposx,...
    pointdwelltime=pointdwelltime,laserpower=laserpowerpf,repetitions=repetitions,dim=2);

sim.defineComponent("estdonut","estimator",@est_qLSQiter2D,parameters={L,probecenter},dim=1:2);
sim.defineComponent("est_x","estimator",@est_quadraticdirect1D,parameters={L},dim=1);
sim.defineComponent("est_y","estimator",@est_quadraticdirect1D,parameters={L},dim=2);

%2D donut
seq={"donut","estdonut"};
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp('no aberration donut:')
sim.summarize_results(out); %display summary of simulation
xs=0:5:65;
displaywhat="rmse";
figure(225);
hold off
sim.scan_fov(seq,xs,clearfigure=true,tag="donut",ax1=displaywhat,linestyle='k',maxlocs=numberOfLocalizations);

%3D donut
seqthnoab={"tophat_xy","estdonut"};
out=sim.runSequence(seqthnoab,"maxlocs",numberOfLocalizations);
disp('no aberration tophat:')
sim.summarize_results(out); %display summary of simulation
sim.scan_fov(seqthnoab,xs, tag="tophat",ax1=displaywhat,linestyle='r',maxlocs=numberOfLocalizations);

% PhaseFLUX
seqpf={"pf_x","est_x","pf_y","est_y"};
out=sim.runSequence(seqpf,"maxlocs",numberOfLocalizations);
disp('no aberration phaseflux:')
sim.summarize_results(out); %display summary of simulation
sim.scan_fov(seqpf,xs, tag="phaseflux",ax1=displaywhat,linestyle='b',maxlocs=numberOfLocalizations);

psf2Ddonut=psf_vec.imagestack("vortex");
psf3Ddonut=psf_vec.imagestack("tophat");
psfPF=psf_vec.imagestack("halfmoonx");
%% Aberrations:
if ~exist("psf_veca","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_veca=PsfVectorial; %simple 2D donut PSF
end
psf_veca.setpinhole("AU",1);
psf_veca.zerooffset=0.000; %true zero
sys_aberr.Zr(1,1)=4;sys_aberr.Zr(1,2)=0;sys_aberr.Zr(1,3)=0.15; %spherical aberrations 
psf_veca.setpar(sys_aberr)
sim.definePattern("donuta", psf_veca, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpowerdonut,repetitions=repetitions)
sim.definePattern("tophat_xya", psf_veca, phasemask="tophat", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpowertophat,repetitions=repetitions,dim=1:2);
sim.definePattern("pf_xa", psf_veca, phasemask="halfmoonx", zeropos=zeroposx,...
    pointdwelltime=pointdwelltime,laserpower=laserpowerpf,repetitions=repetitions,dim=1);
sim.definePattern("pf_ya", psf_veca, phasemask="halfmoony", zeropos=zeroposx,...
    pointdwelltime=pointdwelltime,laserpower=laserpowerpf,repetitions=repetitions,dim=2);

seqdab={"donuta","estdonut"};
out=sim.runSequence(seqdab,"maxlocs",numberOfLocalizations);
% disp('aberration donut:')
% sim.summarize_results(out); %display summary of simulation
sim.scan_fov(seqdab,xs,tag="donut ab",ax1=displaywhat,linestyle='k-.',maxlocs=numberOfLocalizations);

seqthab={"tophat_xya","estdonut"};
out=sim.runSequence(seqthab,"maxlocs",numberOfLocalizations);
disp('aberration tophat:')
% sim.summarize_results(out); %display summary of simulation
sim.scan_fov(seqthab,xs,tag="tophat ab",ax1=displaywhat,linestyle='r-.',maxlocs=numberOfLocalizations);

seqpfa={"pf_xa","est_x","pf_ya","est_y"};
out=sim.runSequence(seqpfa,"maxlocs",numberOfLocalizations);
disp(' aberration phaseflux:')
sim.summarize_results(out); %display summary of simulation
sim.scan_fov(seqpfa,xs, tag="phaseflux ab",ax1=displaywhat,linestyle='b-.',maxlocs=numberOfLocalizations);

psf2Ddonuta=psf_veca.imagestack("vortex");
psf3Ddonuta=psf_veca.imagestack("tophat");
psfPFa=psf_veca.imagestack("halfmoonx");

%% more fluorophores
fc=FlCollection;
fc.add(fl);
fl2=FlStatic;fl2.pos=[0 0 700];
fc.add(fl2);
fl3=FlStatic;fl3.pos=[50 350 900];
fc.add(fl3);
sim.fluorophores=fc;

seq={"donut","estdonut"};
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp('fl donut:')
sim.summarize_results(out); %display summary of simulation

sim.scan_fov(seq,xs,tag="donut fl",ax1=displaywhat,linestyle='k--',maxlocs=numberOfLocalizations);
seqthnoab={"tophat_xy","estdonut"};
out=sim.runSequence(seqthnoab,"maxlocs",numberOfLocalizations);
disp('fl tophat:')
sim.summarize_results(out); %display summary of simulation
sim.scan_fov(seqthnoab,xs, tag="tophat fl",ax1=displaywhat,linestyle='r--',maxlocs=numberOfLocalizations);


out=sim.runSequence(seqpf,"maxlocs",numberOfLocalizations);
disp('fl phaseflux:')
sim.summarize_results(out); %display summary of simulation
sim.scan_fov(seqpf,xs, tag="phaseflux fl",ax1=displaywhat,linestyle='b--',maxlocs=numberOfLocalizations);


%% background
sim.fluorophores=fl;
sim.background=3; %will be multiplied also by laser power factor

seq={"donut","estdonut"};
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
disp('bg donut:')
sim.summarize_results(out); %display summary of simulation
mean(out.bg_photons_gt)

sim.scan_fov(seq,xs,tag="donut bg",ax1=displaywhat,linestyle='k:',maxlocs=numberOfLocalizations);
seqthnoab={"tophat_xy","estdonut"};
out=sim.runSequence(seqthnoab,"maxlocs",numberOfLocalizations);
disp('bg tophat:')
sim.summarize_results(out); %display summary of simulation
mean(out.bg_photons_gt)
sim.scan_fov(seqthnoab,xs, tag="tophat bg",ax1=displaywhat,linestyle='r:',maxlocs=numberOfLocalizations);

out=sim.runSequence(seqpf,"maxlocs",numberOfLocalizations);
disp('bg phaseflux:')
sim.summarize_results(out); %display summary of simulation
sim.scan_fov(seqpf,xs, tag="phaseflux bg",ax1=displaywhat,linestyle='b:',maxlocs=numberOfLocalizations);

ylim([0 L/2])


%% plot PSFs
midp=ceil(size(psf2Ddonut,3)/2);

pa2d=cat(1,psf2Ddonut(:,:,midp),psf2Ddonuta(:,:,midp));
pa3d=cat(1,psf3Ddonut(:,:,midp),psf3Ddonuta(:,:,midp));
papf=cat(1,psfPF(:,:,midp),psfPFa(:,:,midp));

figure(226)
imagesc(cat(2,pa2d/max(pa2d(:)),pa3d/max(pa3d(:)),papf/max(papf(:))))
axis equal