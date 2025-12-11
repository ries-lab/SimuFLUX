function xest=est_octahedron(photonsi,patternposi,L,fwhm,background)
if nargin<5
    background=0;
end
photonsxy=photonsi(1:4);
photonsz=photonsi(5:6);
patternposxy=patternposi(1:4,:);
if length(photonsi)==7 %center probed
    photonsxy(end+1)=photonsi(7);
    photonsz(end+1)=photonsi(7);
    patternposxy(end+1,:)=0;
end
xest=est_donutLSQ1_2D(photonsxy,patternposxy,L,fwhm,background);
xest(3)=est_qLSQiter1D(photonsz,L);

end