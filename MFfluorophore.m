classdef MFfluorophore<handle
    properties
        pos=[0, 0, 0]; %nm
        posmode='diffusion'; %trace, position, function
        brightness=1; %kHz;
        posind=1;
        posparameters=[];
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
        function posout=position(obj,time)
            pos=obj.pos;
            posout=[0 0 0];
            switch obj.posmode
                case 'function'
                    for k=length(pos):-1:1
                        posh(k)=pos{k}(time);
                    end
                case 'position'
                    posh=pos;
                case {'trace','diffusion'}
                    %XXXX add test if ind. already too high (e.g. 2x
                    %position asked. If asked 2x same time: should work)
                    posind=obj.posind;
                    if pos(posind,1)>time
                        posind=find(pos(:,1)<time,1,'last');
                        if isempty(posind)
                            posind=1;
                        end
                    end
                    while pos(posind,1)<time
                        if posind>=size(pos,1)-1 %reached end
                            obj.extendtrace;
                        end
                        posind=posind+1;
                    end
                posh=pos(posind,2:end);
                obj.posind=posind;
            end
            posout(1:length(posh))=posh;
        end
        function makediffusion(obj,D,dt,args)
            arguments
                obj
                D %um^2/s
                dt %us
                args.startpos=[0 0 0]
                args.dim=2;
                args.numpoints=100;
            end
            obj.posparameters={D,dt,args};
            % Dstep=D/1e6*1e6*dt; %D in um^2/s is the same as D in nm^2/us
            Dstep=D*dt;
            time=dt:dt:args.numpoints*dt;
            jumps=(randn(args.numpoints,args.dim)*Dstep);
            pos=horzcat(cumsum(jumps,1))+args.startpos(1:args.dim);
            obj.pos=horzcat(time',pos);
        end
        function extendtrace(obj)
            currentpos=obj.pos(end,:);
            oldpos=obj.pos;
            % obj.posparameters(3).startpos=currentpos(2:end);
            args=obj.posparameters{3};
            obj.makediffusion(obj.posparameters{1},obj.posparameters{2},dim=args.dim,numpoints=args.numpoints,startpos=currentpos(2:end))
            obj.pos(:,1)=obj.pos(:,1)+currentpos(1); %time continues
            obj.pos=vertcat(oldpos,obj.pos);
        end
    end
end
