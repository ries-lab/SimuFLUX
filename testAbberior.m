fl=blinkingfluorophore;
fl.pos=[400,0,0];
fl.brightness=1000;
fl.photonbudget=20000;
psfdonut=PSFMF_donut2D;
sim=AbberiorSimulator(fl);

numlocs=12;
% fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit1000_.json';
fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit100__locLim50_.json';
sim.loadsequence(fname)
% sim.sequence_json.damping=.1;
sim.makepatterns
out=sim.runSequence;
figure(88);hold off;plot(out.loc.loccounter, out.loc.xnm);hold on; plot(out.loc.loccounter,out.loc.xfl);hold on; plot(out.loc.loccounter,out.loc.xgalvo)
plot(out.loc.loccounter,out.loc.xeod)
xlabel('time (itr)')
legend('xest', 'xfl','xgalvo','EOD')