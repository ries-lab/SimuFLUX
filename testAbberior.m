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
% fc.add(fl);
% fl.brightness=1000;
% fl.photonbudget=20000;s
% psfs=PSFMF_donut2D;
sim=AbberiorSimulator(fl);

% psfs=PSFMF_vectorial;
psfpar(1)=struct('psfpar','flat');
psfpar(2)=struct('psfpar','vortex');


numlocs=12;
% fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit1000_.json';
fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit100__locLim50_.json';
sim.loadsequence(fname)
% sim.sequence_json.damping=.1;
sim.makepatterns(psfs,psfpar);
% sim.makepatterns;
% sim.makescoutingpattern([0 0; 1500 1500 ])
out=sim.runSequence;
% out=sim.scoutingSequence(maxrep=1);
figure(88);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl1);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
plot(out.loc.loccounter,out.loc.xeod)
xlabel('time (itr)')
legend('xest', 'xfl','xgalvo','EOD')
sim.displayresults(out,0,50)