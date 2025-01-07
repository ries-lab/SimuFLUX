classdef staticfluorophore<handle
    properties
        pos=[0, 0, 0]; %nm
        brightness=1; %kHz;
    end
    methods
        function Io=intensity(obj,I0,dwelltime, brightness)
            if nargin<4
                brightness=obj.brightness;
            end      
            Io=(brightness)*I0*dwelltime;
        end
        function ph=photons(obj,I0,varargin)
            ph=poissrnd(obj.intensity(I0,varargin{:}));
        end
        function io=tointensity(ii)
            io=ii*obj.brightness+obj.background;
        end
        function reset(obj)
        end
    end
end
