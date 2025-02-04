classdef FlMoving<Fluorophore
    properties
        % pos=[0, 0, 0]; %nm     
        posind=1;
        posparameters=[];
        posmode
    end
    methods
        % function reset(obj)
        %     reset@MFfluorophore(obj);
        %     %later: also reset trace / position if wanted
        % end
        function [posout,isactive]=position(obj,time,props)
            isactive=true;
            pos=obj.pos;
            posout=[0 0 0];
            posmode=obj.posmode;
            switch posmode
                case 'function'
                    for k=length(pos):-1:1
                        posh(k)=pos{k}(time);
                    end
                case 'static'
                    posh=pos;
                case {'trace','diffusion','steps'}
                    %XXXX add test if ind. already too high (e.g. 2x
                    %position asked. If asked 2x same time: should work)
                    posind=max(1,obj.posind-1); 
                    if pos(posind,1)>time
                        posind=find(pos(:,1)<time,1,'last');
                        if isempty(posind)
                            posind=1;
                        end
                    end
                    while pos(posind,1)<time
                        if posind>=size(pos,1)-1 %reached end
                            obj.extendtrace;
                            pos=obj.pos;
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
                args.numpoints=10000;
                args.boundarybox=[]; %nm, half length of box, periodic boundary conditions
            end
            obj.posparameters={D,dt,args};
            % Dstep=D/1e6*1e6*dt; %D in um^2/s is the same as D in nm^2/us.
            % to obtain nm^2/ms multiply by 1000;
            Dstep=D*dt*1000;
            time=dt:dt:args.numpoints*dt;
            jumps=(randn(args.numpoints,args.dim)*Dstep);
            pos=horzcat(cumsum(jumps,1))+args.startpos(1:args.dim);
            if ~isempty(args.boundarybox) %periodic boundary conditions
                if numel(args.boundarybox)==1
                    args.boundarybox=args.boundarybox.*ones(1,args.dim);
                end
                pos2=pos+args.boundarybox;
                pos3=mod(pos2,2*args.boundarybox);
                pos=pos3-args.boundarybox;
            end

            obj.pos=horzcat(time',pos);
            obj.posmode='diffusion';
        end
        function makesteps(obj,stepsize,dwelltime,dt,args)
            arguments
                obj
                stepsize %nm
                dwelltime %ms;
                dt %ms
                args.startpos=[0 0 0]
                args.dim=2;
                args.numpoints=1000;
                args.angle=0;
            end
            obj.posparameters={stepsize,dwelltime,dt,args};
            time=(dt:dt:args.numpoints*dt)';
            
            numjumps=ceil(args.numpoints*dt/(dwelltime)*2+5);
            jmp=exprnd(dwelltime/dt,numjumps,1);
            jumppos=(cumsum(jmp));
            jumppos(jumppos>length(time))=[];
            % xj=zeros(size(time));
            % xj(jumppos)=stepsize;
            xjh=histcounts(jumppos,(0:length(time)))*stepsize;
            xpos=cumsum(xjh');

            R = [cosd(args.angle) -sind(args.angle); sind(args.angle) cosd(args.angle)];
            pos=(R*horzcat(xpos,zeros(size(xpos)))')'+args.startpos(1:args.dim);
            obj.pos=horzcat(time,pos);
            obj.posmode='steps';
        end
        function extendtrace(obj)
            currentpos=obj.pos(end,:);
            oldpos=obj.pos;
            % obj.posparameters(3).startpos=currentpos(2:end);
            args=obj.posparameters{end};
            switch obj.posmode
                case 'diffusion'
                    obj.makediffusion(obj.posparameters{1},obj.posparameters{2},dim=args.dim,numpoints=args.numpoints,startpos=currentpos(2:end),boundarybox=args.boundarybox)
                case 'steps'
                    obj.makesteps(obj.posparameters{1},obj.posparameters{2},obj.posparameters{3},dim=args.dim,numpoints=args.numpoints,startpos=currentpos(2:end), angle=args.angle)
            end        
            obj.pos(:,1)=obj.pos(:,1)+currentpos(1); %time continues
            obj.pos=vertcat(oldpos,obj.pos);
        end
    end
end
