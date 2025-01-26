function [psfs, psfpars]=psf_sequence(sq,psfvec)
% if vectorial in global: use only one PSF object, if also in Itr. other
% parameters are ignored
if isfield(sq.global,'Mode') && strcmp(sq.global.Mode,'PSFMF_vectorial')
    if nargin <2
        psfvec=PSFMF_vectorial;
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
end

for k=1:length(sq.Itr)
    if strcmp(sq.Itr(k).Mode,'PSFMF_vectorial')
        psfs{k}=psfvec;
        psfpars{k}=sq.Itr(k).par;
    else
        psfs{k}=eval(sq.Itr(k).Mode);
        psfpars{k}=sq.Itr(k).par;
    end
end
end