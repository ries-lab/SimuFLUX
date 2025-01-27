%% Vectorial PSF

fl=Fl_static; %define a static fluorophore
fl.pos=[10 0 0];
fl.brightness=1000; %kHz if excited at the center of a Gaussian beam

if ~exist("psf_vec","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec=PSF_vectorial; %simple 2D donut PSF
end
psf_vec.zerooffset=0.000; %true zero

sim=Sim_Simulator(fl); %make a simulator and attach fluorophore

numberOfLocalizations=1000;

%define scan pattern
L=75; %size of scan pattern
orbitpoints=6; %number of probing points in orbit
probecenter=true; %should we also probe the center?
laserpower=5; %relative, increases brightness
pointdwelltime=100; %measurement time in each point
repetitions=2; %how often to repeat the pattern scan
sim.definePattern("donut", psf_vec, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
% sim.definePattern(key, PSF_object, arguments...)

%we need an estimator. Define as component
fwhm=360;% size of the donut, needed for proper estimation. 
sim.defineComponent("estgauss","estimator",@estimators,parameters={"donut2D",sim.patterns("donut").pos,L,fwhm},dim=1:2);
% sim.defineComponent(key,type (estimator),function handle of estimator function,parameters);

%sequence: 
seq={"donut","estgauss"};
out=sim.runSequence(seq);
% out.loc: localizations
% out.fluorophores: position of fluorophores
% out.raw: photon measurements

%calculate CRB for this fluorophore position and scan pattern:
sigmaCRB=sim.calculateCRB("donut",dim=1:2)/sqrt(mean(out.loc.phot));

disp('vectorial PSF:')
sim.displayresults("donut",out); %display summary of simulation

psf0=psf_vec.imagestack("vortex");

%% Zero offset
%now let's add an offset to the PSF to make the minium non-zero
psf_vec.zerooffset=0.005;
out=sim.runSequence(seq);
disp("zero offset = " + psf_vec.zerooffset + ":")
sim.displayresults("donut",out);

%% Aberrations
% let us change the PSF by adding aberrations. Note, in this case we have
% to define the pattern again to calculate the PSFs anew. Instead here, we
% create a second PSF object.

if ~exist("psf_vec2","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vec2=PSF_vectorial; %simple 2D donut PSF
end
%Add Zernike:
% Zr(k,1): n, Zr(k,2): m, Zr(k,3): amplitude as fraction of wavelength
sys.Zr(1,1)=4;sys.Zr(1,2)=0;sys.Zr(1,3)=0.1; %spherical aberrations 
sys.Zr(1,1)=2;sys.Zr(1,2)=2;sys.Zr(1,3)=0.2; %astigmatism 
psf_vec2.setpar(sys)
sim.definePattern("donut_aber", psf_vec2, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
seq={"donut_aber","estgauss"};
out=sim.runSequence(seq);

disp("aberrations:")
sim.displayresults("donut_aber",out);
psfab=psf_vec2.imagestack("vortex");
imx(horzcat(psf0,psfab),'Parent',figure(121)); %compare the two PSFs

%% Pinhole
% We simulate a pinhole in the detection channel
if ~exist("psf_vecph","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vecph=PSF_vectorial; %simple 2D donut PSF
end
psf_vecph.setpinhole("AU",1);
sim.definePattern("donut_ph", psf_vecph, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
seq={"donut_ph","estgauss"};
disp("pinhole:")
out=sim.runSequence(seq);
sim.displayresults("donut_ph",out);

%% Misaligned pinhole
% now lets move the pinhole (misalignment)
if ~exist("psf_vecph2","var") %if PSF is already defined, we need not recalculate it if no parameters are changed
    psf_vecph2=PSF_vectorial; %simple 2D donut PSF
end
psf_vecph2.setpinhole("AU",1,"offset",[150,0]);
sim.definePattern("donut_ph", psf_vecph2, phasemask="vortex", makepattern="orbitscan", orbitpoints=orbitpoints, ...
    probecenter=probecenter,orbitL=L,pointdwelltime=pointdwelltime,laserpower=laserpower,repetitions=repetitions)
seq={"donut_ph","estgauss"};
out=sim.runSequence(seq);
disp("pinhole misaligned:")
sim.displayresults("donut_ph",out);

psfph=psf_vecph2.imagestack("vortex");
imx(horzcat(psf0,psfph),'Parent',figure(123)); %compare the two PSFs