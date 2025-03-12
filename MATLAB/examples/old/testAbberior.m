% clear fl
% xpos=[100 1000];
% % xpos=0;
% for k=1:length(xpos)
%     fl(k)=MFfluorophore;
%     fl(k).pos=[xpos(k),0,0];
% end
% fc=Fluorophorecollection;
fl=MFfluorophore;
fl.pos=[50 0 0];
fl.brightness=10000;
% fc.add(fl);

if ~exist('sim','var') || ~isa(sim,"AbberiorSimulator")
    sim=AbberiorSimulator(fl);
end
sim.posgalvo=[0 0 0];
fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit100__locLim50_.json';
fname2='PSFvectorialJonas.json';
sim.loadsequence(fname,fname2)

numlocs=12;
sim.makepatterns;
% sim.makescoutingpattern([0 0; 1500 1500 ])
out=sim.runSequence;
% out=sim.scoutingSequence(maxrep=1);
figure(88);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl1);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
plot(out.loc.loccounter,out.loc.xeod)
xlabel('time (itr)')
legend('xest', 'xfl','xgalvo','EOD')
sim.displayresults(out,[],50)