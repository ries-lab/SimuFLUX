function out=backgroundsubtractor(out,bg)
bgpoint=out.par.pointdwelltime*bg*out.repetitions;
out.phot=out.phot-bgpoint;
out.bgphot_est=bgpoint;
end