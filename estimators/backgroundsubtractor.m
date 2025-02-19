function out=backgroundsubtractor(out,bg)
bgtot=bg*out.repetitions;
bgpoint=out.pointdwelltime*bgtot;
% out.photrate=out.photrate-bgpoint;
out.photrate=out.photrate-bgtot;
out.bgphot_est=bgpoint;
end