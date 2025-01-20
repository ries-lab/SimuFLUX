classdef PSFMF<handle
    properties
        PSFs=dictionary;
        zerooffset=0; % PSF+zerooffset
        sigmaz=0; %don't do z sectioning.
    end
    methods
        function io=intensity(obj, flpos ,patternpos, phasepattern, L)
       
        end
        function calculatePSFs(obj,phasepattern,Lxs)
           
        end
        function savePSF(obj,name)

        end
        function loadPSF(obj,name)

        end

         function setpinholesimple(obj,args)
            arguments
                obj 
                args.lambda=600; %nm
                args.AU =1; %in airy units
                args.diameter=[] %nm
                args.NA=1.5;
                args.refractiveIndex=1.5; %1.33 water, %1.5 oil
            end
            if isempty(args.diameter)
                args.diameter=args.AU*1.22*args.lambda/args.NA;
            end
            n=args.refractiveIndex;
            fwhm2=(0.88*args.lambda/(n-sqrt(n^2-args.NA^2)))^2+(sqrt(2)*n*args.diameter/args.NA)^2;
            obj.sigmaz=sqrt(fwhm2)/2.35;
         end

         function io=pinholezfac(iin,flposrel)
             sigmaz=obj.sigmaz;
             if sigmaz>0
                 io=exp(-flposrel(:,3).^2/2/sigmaz^2);
             else
                 io=iin;
             end
         end

    end
end

