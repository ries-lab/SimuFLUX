% Run an abberior sequence file and see how EOD an Galvo move.
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist('sim','var') || ~isa(sim,"SimSequencefile")
    sim=SimSequencefile;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];
end
laserpower=1;
fl=FlBleach; %define a bleaching fluorophore
fl.photonbudget=200000;
fl.pos=[200 50 0];
fl.brightness=300; %kHz 
sim.fluorophores=fl;
sim.background=30*0; % a background that is not matched with a proper estimate leads to 'tails'
sim.background_estimated=0; % a similar result is obtained when no background is present but the background is underestimatedd (negative background estimate ). This leads to a bias in the estimator and appearance of "tails"
            %over-estimation of background leads to instabilities

fname='Tracking_2D.json';
fname2='PSFvectorial2D.json'; %use a PSF that is defined via a json file
sim.loadsequence(fname,fname2);

%test different confounding factors, otherwise comment out:
sim.psfvec.setpar('beadradius',0*50e-9) %in m can also lead to tails. Set to zero if no bead used
% sim.sequence.PSF.Itr(1).estimator.par{3}=80; %correct: 120. wrong parameter for first estimator: this makes tails much more pronounced, as currently the Gaussian estimator seems too good.
sim.sequence.locLimit=100; %only track for 1000 localizations
sim.makepatterns;

out=sim.runSequence("repetitions",1);

figure(230);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl1);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
plot(out.loc.loccounter,out.loc.xeod)
xlabel('time (itr)')
ylabel('x position(nm)')
legend('estimated', 'fluorophore','xgalvo','EOD')
indf=out.loc.itr==max(out.loc.itr) & out.loc.vld==1;
sim.summarize_results(out,filter=indf); %display summary of simulation

out.loc.efo(end)