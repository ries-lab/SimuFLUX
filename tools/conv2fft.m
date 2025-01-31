function c=conv2fft(a,b)
afft=fft2(a);
bfft=fft2(b,size(a,1),size(a,2));
cfft=afft.*bfft;
c=real(ifft2(cfft));
for k=1:2 %dims
    c=ifftshift(c,k);
end
