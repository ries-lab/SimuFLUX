classdef bleachingFluorophore<MFfluorophore
    properties
        photonbudget=inf; %1000 photons
    end
    methods
        function Io=intensity(obj,I0,dwelltime,phfac)
            intensity=(obj.brightness/1000)*I0*dwelltime;
            remainingphotons=obj.remainingphotons;
            Io=min(intensity,remainingphotons);
            obj.remainingphotons=remainingphotons-Io/phfac;
        end
        
        function reset(obj) %new fluorophore
           
            % recalculate photon budget
            obj.remainingphotons=exprnd(obj.photonbudget);
        end
 
    end
end

