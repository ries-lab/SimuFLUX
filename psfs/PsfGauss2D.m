classdef PsfGauss2D<Psf
    properties
        % fwhm=310; %comaprison with calculated PSF
        sigma=310/2.35; %310 from donut
    end
    methods
        function [io,phfac]=intensity(obj, flpos ,patternpos, phasepattern, L)
            flposrel=flpos-patternpos;
            r2=flposrel(:,1:2).^2;
            sigma=obj.sigma;
            rs2=sum(r2,2);
         
            % zerooffset=obj.zerooffset;
            io=exp(-rs2/2/sigma^2);%+zerooffset; zerooffset not important for Gauss
            if obj.sigmaz>0
                phfac=obj.pinholezfac(flposrel);
                io=io.*phfac;
            else
                phfac=ones(size(io));
            end
            
        end
        

        function plotprofile(obj)
            xaxv=-700:700;
            pos(:,1)=xaxv;
            pos(:,3)=0;
            inten=obj.intensity(pos,0,0);
            figure(33)
            plot(xaxv, 0.3*inten(:,round(end/2),round(end/2)))
        end
    end
end
