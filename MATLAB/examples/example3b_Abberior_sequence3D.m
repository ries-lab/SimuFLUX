% Run an abberior sequence file and see how EOD an Galvo move.
addpath(genpath(fileparts(fileparts(mfilename('fullpath'))))); %add all folders to serach path
if ~exist('sim','var') || ~isa(sim,"SimSequencefileAbberior")
    sim=SimSequencefileAbberior;
else
    sim.posgalvo=[0 0 0];sim.posEOD=[0 0 0];
end
laserpower=1;
fl=FlBleach; %define a bleaching fluorophore
fl.photonbudget=20000000;
fl.pos=[140 10 100];
fl.brightness=3000; %kHz 
sim.fluorophores=fl;
sim.background=30*0; % a background that is not matched with a proper estimate leads to 'tails'
sim.background_estimated=0; % a similar result is obtained when no background is present but the background is underestimatedd (negative background estimate ). This leads to a bias in the estimator and appearance of "tails"
            %over-estimation of background leads to instabilities
sim.psfvec.setpinhole("AU",1) %Abberior: imspector pinhole overwrites the pinhole in the settings file

fname='Tracking_3D.json';
% fname='Imaging_3D.json';

sim.loadsequence(fname);

%test different confounding factors, otherwise comment out:
sim.psfvec.setpar('beadradius',0*50e-9) %in m can also lead to tails. Set to zero if no bead used
sim.sequence.locLimit=100; %only track for 1000 localizations
sim.makepatterns;

out=sim.runSequence("repetitions",1);

figure(230);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl1);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
plot(out.loc.loccounter,out.loc.xeod)
xlabel('time (itr)')
ylabel('x position(nm)')
legend('estimated', 'fluorophore','xgalvo','EOD')

figure(231);hold off;plot(out.loc.loccounter, out.loc.znm);hold on; plot(out.loc.loccounter,out.loc.zfl1);hold on; plot(out.loc.loccounter,out.loc.zgalvo)
plot(out.loc.loccounter,out.loc.zeod)
xlabel('time (itr)')
ylabel('x position(nm)')
legend('estimated', 'fluorophore','xgalvo','EOD')


indf=out.loc.itr>=max(out.loc.itr)+sim.sequence.headstart+1 & out.loc.vld==1;
sim.summarize_results(out,filter=indf); %display summary of simulation

out.loc.efo(end)