classdef staticfluorophore<handle
    properties
        pos=[0, 0, 0]; %nm
        brightness=1; %kHz;
        blinking=false;
        ton
        toff
    end
    methods
        function Io=intensity(obj,I0,brightness)
            if nargin<3
                brightness=obj.brightness;
            end      
            Io=(brightness)*I0;
        end
        function ph=photons(obj,I0,varargin)
            ph=poissrnd(obj.intensity(I0,varargin{:}));
        end
        function io=tointensity(ii)
            io=ii*obj.brightness+obj.background;
        end
    end
end
