classdef MFfluorophore<handle
    properties
        pos=[0, 0, 0]; %nm
        brightness=1000; %kHz;
        remainingphotons=inf;
        %dummy properties for compatibility with collection
        numberOfFluorophores=1;
        allbleached=false;
    end
    methods
        function Io=intensity(obj,I0,dwelltime, phfac)
            %dwelltime: us, brightness kHz
            % if nargin<4
            %     brightness=obj.brightness;
            % end      
            Io=(obj.brightness/1000)*I0*dwelltime.*phfac;
        end
        function ph=photons(obj,I0,varargin)
            ph=poissrnd(obj.intensity(I0,varargin{:}));
        end
        % function io=tointensity(obj,ii)
        %     io=ii*(obj.brightness/10000);
        % end
        function reset(obj)
        end
        function [posout,isactive]=position(obj,time)
            posout=obj.pos;
            isactive=true;
        end

        %dummy functions for compatibility with collection
        function updateonoff(obj,time)
        end
    end
end
