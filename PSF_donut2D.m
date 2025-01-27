classdef PSF_donut2D<PSF_Pointspreadfunction
    properties
        % fwhm=310; %comaprison with calculated PSF
    end
    methods
        function [io,phfac]=intensity(obj, flpos ,patternpos, phasepattern, L)
            flposrel=flpos-patternpos;
            r2=flposrel(:,1:2).^2;
            fwhm=obj.fwhm;
            rs2=sum(r2,2)/fwhm^2;
            zerooffset=obj.zerooffset;
            % ioa=0.3*4*exp(1)*log(2)*rs2.*exp(-4*log(2)*rs2)+zerooffset;    %0.3: comparison with calculated PSF 
            io=2.2610*rs2.*exp(-2.7726*rs2)+zerooffset;
            if obj.sigmaz>0
                phfac=obj.pinholezfac(flposrel);
                io=io*phfac;
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
