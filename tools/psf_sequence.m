function [psfs, phasemasks]=psf_sequence(sq,psfvec,seq)
% if vectorial in global: use only one PSF object, if also in Itr. other
% parameters are ignored
if nargin>2 %abberior iteration: extrac pinhole from iterations, overwrite default.
    sq.global.addpar.pinhole={'AU',seq.Itr(end).Mode.phDiaAU};
end
if isfield(sq,'global') && isfield(sq.global,'Mode') && strcmp(sq.global.Mode,'PSF_vectorial')
    if nargin <2
        psfvec=PSF_vectorial;
    end
    fng=fieldnames(sq.global);
    % fng=intersect(fng,{'sys','opt','addpar'});
    
    for k=1:length(fng)
        ssth=sq.global.(fng{k});
        if isstruct(ssth)
            psfvec.setpar(sq.global.(fng{k}));
        else
            psfvec.setpar(struct(fng{k},sq.global.(fng{k})))
        end
           
    end
    if isfield(sq.global.addpar,"pinhole")
        psfvec.setpinhole(sq.global.addpar.pinhole{:})
    end
end

for k=1:length(sq.Itr)
    if strcmp(sq.Itr(k).Mode,'PSF_vectorial')
        psfs{k}=psfvec;
        phasemasks{k}=sq.Itr(k).par;
    else
        psfs{k}=eval(sq.Itr(k).Mode);
        phasemasks{k}=sq.Itr(k).par;
    end
end
end