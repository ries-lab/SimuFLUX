classdef MFfluorophore<handle
    properties
        pos=[0, 0, 0]; %nm
        brightness=1000; %kHz;
        remainingphotons=inf;
    end
    methods
        function Io=intensity(obj,I0,dwelltime, phfac)
            %dwelltime: us, brightness kHz
            % if nargin<4
            %     brightness=obj.brightness;
            % end      
            Io=(obj.brightness/1000)*I0*dwelltime*phfac;
        end
        function ph=photons(obj,I0,varargin)
            ph=poissrnd(obj.intensity(I0,varargin{:}));
        end
        function io=tointensity(ii)
            io=ii*(obj.brightness/10000)+obj.background;
        end
        function reset(obj)
        end
        function posout=position(obj,time)
            posout=obj.pos;
        end
    end
end
