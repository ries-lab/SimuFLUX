classdef Fl_bleach<Fl_Fluorophore
    properties
        photonbudget=inf; %1000 photons
    end
    methods
        function Io=intensity(obj,I0,dwelltime,time,phfac,props)
            remainingphotons=obj.remainingphotons;
            if remainingphotons==0
                Io=0;
                return
            end
            intensity=(obj.brightness/1000)*I0*dwelltime;
            
            Io=min(intensity,remainingphotons);
            obj.remainingphotons=remainingphotons-Io/phfac;
        end
        function set.photonbudget(obj,val)
            obj.photonbudget=val;
            obj.reset(0);
        end
        
        function reset(obj,time) %resets the counter for bleaching, i.e. switches the fluorophore on again
            obj.remainingphotons=exprnd(obj.photonbudget);
        end
    end
end

