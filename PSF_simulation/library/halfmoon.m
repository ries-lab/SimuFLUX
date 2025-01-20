%% This function introduces a phase shift at half the area of the input image.
% Code written by Takahiro DEGUCHI, European Molecular Biology Laboratory, October 2023.
% This is meant to be used with MATLAB library, "Electromagnetic field in the focus region of a microscope objective", 
% based on the original publication, Leutenegger et al., "Fast focus field calculations", published in Optics Express 14 (2006) 


function E = halfmoon(sys,E,r,t,p)
[xx,yy] = pol2cart(p,r);
xx=xx+sys.maskshift(1);
yy=yy+sys.maskshift(2);
[p,r]=cart2pol(xx,yy);

idx = p>=0;       
idx = 2*idx-1;
phase = idx.*exp((pi)*i);   
phase(p>=0) = phase(p>=0).*exp(sys.delshift*i);
phase = repmat(phase,size(E,1),1);
E = E.*phase;
end

