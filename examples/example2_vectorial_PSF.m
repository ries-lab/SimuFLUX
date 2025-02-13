%% Vectorial PSF
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path

fl=FlStatic(brightness=1000); %define a static fluorophore
fl.pos=[10 0 0];

if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PsfVectorial; %simple 2D donut PSF
end
psf_vec.zerooffset=0.000; %true zero

sim=Simulator(fluorophores=fl); %make a simulator and attach fluorophore

numberOfLocalizations=1000;

%define scan pattern
L=75; %size of scan pattern
probecenter=true; %should we also probe the center?
orbitpoints=6;
laserpower=5; %relative, increases brightness
pointdwelltime=0.1; % ms, measurement time in each point
repetitions=2; %how often to repeat the pattern scan
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)

%we need an estimator. Define as component
sim.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
% sim.defineComponent(key,type (estimator),function handle of estimator function,parameters);

%sequence: 
seq={"donut","estdonut"};
out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
% out.loc: localizations
% out.fluorophores: position of fluorophores
% out.raw: photon measurements

%calculate CRB for this fluorophore position and scan pattern:
sigmaCRB=sim.calculateCRB("donut",dim=1:2)/sqrt(mean(out.loc.phot));

disp('vectorial PSF:')
sim.summarize_results(out); %display summary of simulation

psf0=psf_vec.imagestack("vortex");

%% Zero offset
%now let's add an offset to the PSF to make the minium non-zero
psf_vec.zerooffset=0.005;
out=sim.runSequence(seq);
disp("zero offset = " + psf_vec.zerooffset + ":")
sim.summarize_results(out);

%% Aberrations
% let us change the PSF by adding aberrations. Note, in this case we have
% to define the pattern again to calculate the PSFs anew. Instead here, we
% create a second PSF object.

if ~exist("psf_vec2","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec2=PsfVectorial; %simple 2D donut PSF
end
%Add Zernike:
% Zr(k,1): n, Zr(k,2): m, Zr(k,3): amplitude as fraction of wavelength
sys_aberr.Zr(1,1)=4;sys_aberr.Zr(1,2)=0;sys_aberr.Zr(1,3)=0.1; %spherical aberrations 
sys_aberr.Zr(2,1)=2;sys_aberr.Zr(2,2)=2;sys_aberr.Zr(2,3)=0.05; %astigmatism 
psf_vec2.setpar(sys_aberr)
sim.definePattern("donut_aber", psf_vec2, phasemask="vortex", makepattern="orbitscan", orbitpoints=4, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
seq={"donut_aber","estdonut"};
out=sim.runSequence(seq);

disp("aberrations:")
sim.summarize_results(out);
psfab=psf_vec2.imagestack("vortex");
imx(horzcat(psf0,psfab),'Parent',figure(220),'Title',"Aberrations"); %compare the two PSFs

%% Misaligned phase plate
if ~exist("psf_vec2","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec2=PsfVectorial; %simple 2D donut PSF
end
%Add Zernike:
% Zr(k,1): n, Zr(k,2): m, Zr(k,3): amplitude as fraction of wavelength
sys_mis.Zr(1,1)=4;sys_mis.Zr(1,2)=0;sys_mis.Zr(1,3)=0.0; %spherical aberrations 
sys_mis.Zr(2,1)=2;sys_mis.Zr(2,2)=2;sys_mis.Zr(2,3)=0.0; %astigmatism 
sys_mis.maskshift=[0.2,0]; % radius of pupil function is 1
psf_vec2.setpar(sys_mis)
sim.definePattern("donut_misaligned", psf_vec2, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
seq={"donut_misaligned","estdonut"};
out=sim.runSequence(seq);

disp("misaligned phase plate:")
sim.summarize_results(out);
psfab=psf_vec2.imagestack("vortex");
imx(horzcat(psf0,psfab),'Parent',figure(221),'Title',"Misaligned phase plate"); %compare the two PSFs


%% Pinhole
% We simulate a pinhole in the detection channel
if ~exist("psf_vecph","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vecph=PsfVectorial; %simple 2D donut PSF
end
psf_vecph.setpinhole("AU",1);
sim.definePattern("donut_ph", psf_vecph, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
seq={"donut_ph","estdonut"};
disp("pinhole:")
out=sim.runSequence(seq);
sim.summarize_results(out);

%% Misaligned pinhole
% now lets move the pinhole (misalignment)
if ~exist("psf_vecph2","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vecph2=PsfVectorial; %simple 2D donut PSF
end
psf_vecph2.setpinhole("AU",1,"offset",[150,0]);
sim.definePattern("donut_ph", psf_vecph2, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
seq={"donut_ph","estdonut"};
out=sim.runSequence(seq);
disp("pinhole misaligned:")
sim.summarize_results(out);

psfph=psf_vecph2.imagestack("vortex");
imx(horzcat(psf0,psfph),'Parent',figure(222),'Title',"Misaligned pinhole"); %compare the two PSFs


%% 3D with tophat
sim.fluorophores.pos=[0,0, 50];
if ~exist("psf_vecth","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vecth=PsfVectorial; %simple 2D donut PSF
end
psf_vecth.setpinhole("AU",1);
sys_th.Zr(1,1)=4;sys_th.Zr(1,2)=0;sys_th.Zr(1,3)=0.0; %spherical aberrations 
sys_th.Zr(2,1)=2;sys_th.Zr(2,2)=2;sys_th.Zr(2,3)=0.0; %astigmatism 
psf_vecth.setpar(sys_th)

orbitpoints=4;
probecenterxy=true;
probecenterz=true;
L=75;
Lz=150;
% sigmaz=200;
laserpower=30;

sim.definePattern("tophat_xy", psf_vecth, phasemask="tophat", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenterxy,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions,dim=1:2);
sim.definePattern("tophat_z", psf_vecth, phasemask="tophat", makepattern="zscan", orbitpoints=2, ...
    probecenter=probecenterz,orbitL=Lz,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions,dim=3);

sim.defineComponent("esttophat_xy","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
sim.defineComponent("esttophat_z","estimator",@est_quad1Diter,parameters={Lz},dim=3);

seq={"tophat_xy","esttophat_xy","tophat_z","esttophat_z"};
out=sim.runSequence(seq);

disp("3D with tophat:")
sim.summarize_results(out);

% 3D with tophat and vortex
sim.definePattern("donut_xy", psf_vecth, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)

seq={"donut_xy","estdonut","tophat_z","esttophat_z"};
out=sim.runSequence(seq);

disp("3D with donut and tophat:")
sim.summarize_results(out);


%% PhaseFlux 3D localization
sim.fluorophores.pos=[20,-10, 50];
if ~exist("psf_vecphaseflux","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vecphaseflux=PsfVectorial; %simple 2D donut PSF
end
psf_vecphaseflux.setpinhole("AU",1);


L=75;
Lz=150;
fwhm=450;
sigmaz=200;
laserpower=30;
zeroposx=[-1;1;0]*L/2;
zeroposz=[-1;1;0]*Lz/2;

sim.definePattern("pf_x", psf_vecphaseflux, phasemask="halfmoonx", zeropos=zeroposx,...
    pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions,dim=1);
sim.definePattern("pf_y", psf_vecphaseflux, phasemask="halfmoony", zeropos=zeroposx,...
    pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions,dim=2);
sim.definePattern("pf_z", psf_vecphaseflux, phasemask="tophat", zeropos=zeroposz,...
    pointdwelltime=pointdwelltime,laserpower=laserpower*10,repetitions=repetitions,dim=3);

sim.defineComponent("est_x","estimator",@est_quad1Diter,parameters={L},dim=1);
sim.defineComponent("est_y","estimator",@est_quad1Diter,parameters={L},dim=2);
sim.defineComponent("est_z","estimator",@est_quad1Diter,parameters={Lz},dim=3);
% sim.defineComponent("esttophat_z","estimator",@est_quad1Diter,parameters={Lz},dim=3);


seq={"pf_x","est_x","pf_y","est_y","pf_z","est_z"};
% seq={"pf_z","est_z"};
out=sim.runSequence(seq);
disp("PhaseFLUX:")
sim.summarize_results(out);
