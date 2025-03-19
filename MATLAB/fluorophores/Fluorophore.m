classdef Fluorophore<handle
    properties
        pos=[0, 0, 0]; %nm
        brightness=1000; %kHz;
        remainingphotons=inf;
        %dummy properties for compatibility with collection
        numberOfFluorophores=1;
        allbleached=false;
    end
    methods
        function obj=Fluorophore(args)
            arguments
                args.pos=[0,0,0];
                args.brightness=1000;
            end
            obj.pos=args.pos;
            obj.brightness=args.brightness;
        end
        function Io=intensity(obj,I0,dwelltime, time,phfac)
            % if nargin>=6 && ~isempty(props) %for performance: avoid property access
                % brightness=props.brightness;
            % else
                brightness=obj.brightness;
            % end
            
            Io=(brightness)*I0*dwelltime;%.*phfac;
        end
        function ph=photons(obj,I0,varargin)
            ph=poissrnd(obj.intensity(I0,varargin{:}));
        end
        function reset(obj,time)
        end
        function [posout,isactive]=position(obj,time)
            pos=obj.pos;
            posout=pos;
            isactive=true;
        end

        %dummy functions for compatibility with collection
        function updateonoff(obj,time)
        end
        function props=getproperties(obj)
            props.brightness=obj.brightness;
            props.pos=obj.pos;
        end
    end
end
