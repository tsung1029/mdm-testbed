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
read = 1;
fourier = 1;
plot_real = 0;
plot_four = 0;
e_plot = 0;
b_plot = 0;
s_plot = 0;
read_start = 1;
read_stop = 1;
snap = 1;

disp('setting up time and space tables')
%set up temporal tables/ data.
steps = round(tend/(nt*dt));
temporal = [0:steps - 1];
time = nt*dt*temporal;
dw = 2*pi/tend;
omega = dw*(temporal - steps/2);

%steps2 = round((t_stop-t_start)/nt);
%temporal2 = [0:steps2 - 1];
%dw2 = 2*pi/(dt*t_stop-dt*t_start);
%omega2 = dw2*(temporal2 - steps2/2);

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

affp = 2^(indx + indy)/np;

if read == 1

  disp('read in files from 2-D UPIC')

  %read in field data from bbeps2
  disp(' -field data')
  %ex_txy = zeros(2^indx,2^indy,steps);
  %ey_txy = zeros(2^indx,2^indy,steps);
%  ez_txy = zeros(2^indx,2^indy,steps);
%  bx_txy = zeros(2^indx,2^indy,steps);
%  by_txy = zeros(2^indx,2^indy,steps);
  %bz_txy = zeros(2^indx,2^indy,steps);
%  energy = zeros(steps);

  for i = read_start:nt:read_stop

    step = i-1

%    ex_txy(:,:,step+1) = emma_2d('MS','e1',1,step,indx,indy,2);
%    ey_txy(:,:,step+1) = emma_2d('MS','e2',1,step,indx,indy,2);
    ez_txy(:,:,step+1) = emma_2d('MS','e3',1,step,indx,indy,2);

    %bx_txy(:,:,step+1) = emma_2d('MS','b1',1,step,indx,indy,2);
%    by_txy(:,:,step+1) = emma_2d('MS','b2',1,step,indx,indy,2);
%    bz_txy(:,:,step+1) = emma_2d('MS','b3',1,step,indx,indy,2);

  end

  %ex_txy = double(ex_txy);
  %ey_txy = double(ey_txy);
%  ez_txy = double(ez_txy);

  %bx_txy = double(bx_txy);
  %by_txy = double(by_txy);
  %bz_txy = double(bz_txy);

end

if snap > 0

 % ex_xy = squeeze(ex_txy(:,:,snap));
 % ey_xy = squeeze(ey_txy(:,:,snap));
  ez_xy = squeeze(ez_txy(:,:,snap));
 % emag_xy = sqrt(ex_xy.^2 + ey_xy.^2);

 % bx_xy = squeeze(bx_txy(:,:,snap));
 % by_xy = squeeze(by_txy(:,:,snap));
 % bz_xy = squeeze(bz_txy(:,:,snap));

 % sx_xy = ey_xy.*bz_xy - by_xy.*ez_xy;
 % sy_xy = ez_xy.*bx_xy - bz_xy.*ex_xy;
 % sz_xy = ex_xy.*by_xy - bx_xy.*ey_xy;
 % smag = sqrt(sx_xy.^2 + sy_xy.^2);

end

if fourier == 1 && snap > 0

 % ex_kxky = fftshift(fft(fftshift(fft(ex_xy,[],1),1),[],2),2);
 % ey_kxky = fftshift(fft(fftshift(fft(ey_xy,[],1),1),[],2),2);
  ez_kxky = fftshift(fft(fftshift(fft(ez_xy,[],1),1),[],2),2);
 % emag_kxky = sqrt(abs(ex_kxky).^2 + abs(ey_kxky).^2);

 % bx_kxky = fftshift(fft(fftshift(fft(bx_xy,[],1),1),[],2),2);
 % by_kxky = fftshift(fft(fftshift(fft(by_xy,[],1),1),[],2),2);
 % bz_kxky = fftshift(fft(fftshift(fft(bz_xy,[],1),1),[],2),2);
 % bmag_kxky = sqrt(abs(bx_kxky).^2 + abs(by_kxky).^2);

 % sx_kxky = ey_kxky.*bz_kxky - by_kxky.*ez_kxky;
 % sy_kxky = ez_kxky.*bx_kxky - bz_kxky.*ex_kxky;
 % sz_kxky = ex_kxky.*by_kxky - bx_kxky.*ey_kxky;
 % smagk = sqrt(abs(sx_kxky).^2 + abs(sy_kxky).^2);

end



if snap > 0 && plot_real == 1
  disp('  -snapshot in real space')

  [posx,posy] = meshgrid(x,y);

  if e_plot == 1
%    figure
%    pcolor(posx,posy,log(abs(ex_xy)).'); shading flat
%    title(strjoin({'E_x(x,y,t=',num2str(snap),')'}))

%    figure
%    pcolor(posx,posy,log(abs(ey_xy)).'); shading flat
%    title(strjoin({'E_y(x,y,t=',num2str(snap),')'}))

    figure
    %pcolor(posx,posy,log(abs(ez_xy)).'); shading flat
    imagesc(log(abs(ez_xy)).'); shading flat
    title(strjoin({'E_z(x,y,t=',num2str(snap),')'}))

%    figure
%    pcolor(posx,posy,log(emag_xy).'); shading flat
%    title(strjoin({'|E(x,y,t=',num2str(snap),')|'}))

  end

  if s_plot == 1

    figure
    pcolor(posx,posy,log(abs(sx_xy))); shading flat
    title(strjoin({'S_x(x,y,t=',num2str(snap),')'}))

    figure
    pcolor(posx,posy,log(abs(sy_xy))); shading flat
    title(strjoin({'S_y(x,y,t=',num2str(snap),')'}))

    figure
    pcolor(posx,posy,log(smag)); shading flat
    title(strjoin({'|S(x,y,t=',num2str(snap),')|'}))

  end

  if b_plot == 1

    figure
    pcolor(posx,posy,log(abs(bx_xy))); shading flat
    title(strjoin({'B_x(x,y,t=',num2str(snap),')'}))

    figure
    pcolor(posx,posy,log(abs(by_xy))); shading flat
    title(strjoin({'B_y(x,y,t=',num2str(snap),')'}))

%    figure
%    pcolor(posx,posy,log(abs(bz_xy)).'); shading flat
%    title(strjoin({'B_z(x,y,t=',num2str(snap),')'}))

  end

end

if snap > 0 && plot_four == 1
  disp('  -snapshot in Fourier space')

  [poskx,posky] = meshgrid(kx,ky);

  if e_plot == 1
    figure
    pcolor(poskx,posky,log(abs(ex_kxky)).'); shading flat
    title(strjoin({'E_x(kx,ky,t=',num2str(snap),')'}))

    figure
    pcolor(poskx,posky,log(abs(ey_kxky)).'); shading flat
    title(strjoin({'E_y(kx,ky,t=',num2str(snap),')'}))

    figure
    pcolor(poskx,posky,log(abs(ez_kxky)).'); shading flat
    title(strjoin({'E_z(kx,ky,t=',num2str(snap),')'}))

%    figure
%    pcolor(poskx,posky,log(emag_kxky).'); shading flat
%    title(strjoin({'|E(kx,ky,t=',num2str(snap),')|'}))

  end

  if b_plot == 1

    figure
    pcolor(poskx,posky,log(abs(bx_kxky)).'); shading flat
    title(strjoin({'B_x(kx,ky,t=',num2str(snap),')'}))

    figure
    pcolor(poskx,posky,log(abs(by_kxky)).'); shading flat
    title(strjoin({'B_y(kx,ky,t=',num2str(snap),')'}))

    figure
    pcolor(poskx,posky,log(abs(bz_kxky)).'); shading flat
    title(strjoin({'B_z(kx,ky,t=',num2str(snap),')'}))

    figure
    pcolor(poskx,posky,log(bmag_kxky).'); shading flat
    title(strjoin({'|B(kx,ky,t=',num2str(snap),')|'}))

  end

  if s_plot == 1

    figure
    pcolor(poskx,posky,log(abs(sx_kxky))); shading flat
    title(strjoin({'S_x(kx,ky,t=',num2str(snap),')'}))

    figure
    pcolor(poskx,posky,log(abs(sy_kxky))); shading flat
    title(strjoin({'S_y(kx,ky,t=',num2str(snap),')'}))

    figure
    pcolor(poskx,posky,log(smagk)); shading flat
    title(strjoin({'|S(kx,ky,t=',num2str(snap),')|'}))

  end

end

