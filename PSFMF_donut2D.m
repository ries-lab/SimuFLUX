classdef PSFMF_donut2D<PSFMF
    properties
        fwhm=310; %comaprison with calculated PSF
    end
    methods
        function io=intensity(obj,phasepattern, L, flposrel)
            r2=flposrel(:,1:2).^2;
            rs2=sum(r2,2)/obj.fwhm^2;
            io=0.3*4*exp(1)*log(2)*rs2.*exp(-4*log(2)*rs2)+obj.zerooffset;    %0.3: comparison with calculated PSF 
        end
        function calculatePSFs(obj,phasepattern,Lxs)
        end
        function savePSF(obj,name)
        end
        function loadPSF(obj,name)
        end

        function plotprofile(obj)
            xaxv=-700:700;
            pos(:,1)=xaxv;
            pos(:,3)=0;
            inten=obj.intensity(0,0,pos);
            figure(33)
            plot(xaxv, 0.3*inten(:,round(end/2),round(end/2)))
        end
    end
end
