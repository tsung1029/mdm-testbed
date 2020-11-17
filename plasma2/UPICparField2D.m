%clear all

%Basic plotting diagnostic for 2D UPIC
%parameters from 2D UPIC input file
indx = 9;
indy = 9;
tend = 10.5;
dt = 0.1;
nt = 20;
np = 1;
dx = 1;
dy = 1;

%options

fourier = 1;
snap = 1;

disp('setting up space tables')

%spatial tables/data
modesxe = 2^(indx-1);
dkx = 2*pi/(2^indx);
kx = dkx*(-modesxe:modesxe-1);
dx = 2^indx/(2*modesxe);
x = dx*(0:2*modesxe-1);

modesye = 2^(indy-1);
if modesye == 1
  dky = 0.0;
else
  dky = 2*pi/(2^indy);
end
ky = dky*(-modesye:modesye-1);
dy = 2^indy/(2*modesye);
y = dy*(0:2*modesye-1);

if snap > 0

 % ex_xy = squeeze(ex_txy(:,:,snap));
 % ey_xy = squeeze(ey_txy(:,:,snap));
  ez_xy = squeeze(ez_txy(:,:,snap));
 % emag_xy = sqrt(ex_xy.^2 + ey_xy.^2);

 % bx_xy = squeeze(bx_txy(:,:,snap));
 % by_xy = squeeze(by_txy(:,:,snap));
 % bz_xy = squeeze(bz_txy(:,:,snap));

end

if fourier == 1 && snap > 0

 % ex_kxky = fftshift(fft(fftshift(fft(ex_xy,[],1),1),[],2),2);
 % ey_kxky = fftshift(fft(fftshift(fft(ey_xy,[],1),1),[],2),2);
  ez_kxky = fftshift(fft2(ez_xy));
 % emag_kxky = sqrt(abs(ex_kxky).^2 + abs(ey_kxky).^2);

 % bx_kxky = fftshift(fft(fftshift(fft(bx_xy,[],1),1),[],2),2);
 % by_kxky = fftshift(fft(fftshift(fft(by_xy,[],1),1),[],2),2);
 % bz_kxky = fftshift(fft(fftshift(fft(bz_xy,[],1),1),[],2),2);
 % bmag_kxky = sqrt(abs(bx_kxky).^2 + abs(by_kxky).^2);

end

% get the form factor to get E parallel
for ix = 1 : 2*modesxe
  for jy = 1 : 2*modesye
  
    ffc_x(ix,jy) = abs(kx(ix))/sqrt(kx(ix)^2 + ky(jy)^2);
    ffc_y(ix,jy) = abs(ky(jy))/sqrt(kx(ix)^2 + ky(jy)^2);
    
  end
end

ffc_x(modesxe+1,modesye+1) = 1.0;
ffc_y(modesxe+1,modesye+1) = 1.0;
ffc_xp = ffc_x;
ffc_xp(1:modesxe,:) = 0.0;
ffc_xm(modesxe+1:2*modesxe,:) = 0.0;

%exx_kxky = ex_kxky.*ffc_x;
%exy_kxky = ex_kxky.*ffc_y;
ezxp_kxky = ez_kxky.*ffc_xp;
ezx_kxky = ez_kxky.*ffc_x;
ezy_kxky = ez_kxky.*ffc_y;

% take real part because of small imaginary part from roundoff
%exx_xy = fft2(fftshift(exx_kxky));
%exy_xy = fft2(fftshift(exy_kxky));
ezxp_xy = real(fft2(fftshift(ezxp_kxky)))/2^(indx+indy);
ezx_xy = real(fft2(fftshift(ezx_kxky)))/2^(indx+indy);
ezy_xy = real(fft2(fftshift(ezy_kxky)))/2^(indx+indy);



%exx_xy = fftshift(fft(fftshift(fft(exx_kxky,[],1),1),[],2),2);
%exy_xy = fftshift(fft(fftshift(fft(exy_kxky,[],1),1),[],2),2);


      

