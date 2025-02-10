%% Vectorial PSF
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path

fl=FlStatic(brightness=1000); %define a static fluorophore
fl.pos=[0 0 0];



pixelsizes=[5 10 15 20 30 50 100 150]* 1e-9;
pixelsizes=pixelsizes(end:-1:1);
ix2nm=load("examples/donut2nmpixelprofile.mat").ix;
xs=0:0.01:20;
figure(88)
hold off
pixelsizes=150e-9;
for k=1:length(pixelsizes)
    psf_vec=PsfVectorial; %simple 2D donut PSF
    psf_vec.setpar('dr',pixelsizes(k),'dz',pixelsizes(k))
    
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
    sim.defineComponent("estdonut","estimator",@est_quad2Diter,parameters={L,probecenter},dim=1:2);
    seq={"donut","estdonut"};
    out=sim.runSequence(seq,"maxlocs",numberOfLocalizations);
    
    disp('vectorial PSF:')
    sim.summarize_results(out); %display summary of simulation
    
    psf0=psf_vec.imagestack("vortex");
   
    pos(:,1)=xs;
    pos(:,3)=0;
    ix=psf_vec.intensity(pos,[0 0 0],"vortex0");
    loglog(xs,(ix))
    hold on
end
loglog(xs,ix2nm)
legend([string(pixelsizes) string(2)])
xlabel('x position (nm)')
ylabel('Intensity)')
title("sampling error PSF")
% save("examples/donut2nmpixelprofile.mat","ix");