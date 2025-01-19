fl=blinkingfluorophore;
fl.brightness=1000;
fl.photonbudget=100000;
psfdonut=PSFMF_donut2D;
sim=AbberiorSimulator(fl);
numlocs=12;
% fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit1000_.json';
fname='/Users/ries/Downloads/2Dtracking_L75nm_phtLimit100__locLim50_.json';
sim.loadsequence(fname)
sim.makepatterns
out=sim.runSequence;