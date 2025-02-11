function out=backgroundsubtractor(out,bg)
bgpoint=out.pointdwelltime*bg*out.repetitions;
out.photrate=out.photrate-bgpoint;
out.bgphot_est=bgpoint;
end