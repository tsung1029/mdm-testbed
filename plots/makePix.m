read_start = 0;
read_stop = 500;
nt = 10;

for i = read_start:nt:read_stop

%    frame_index = i-1
     frame_index = i

scrsize = get(0,'Screensize');
      figure('Visible','off');set(gcf, 'Position', [scrsize(1) scrsize(2) scrsize(3) floor(scrsize(4))]);
      set(gcf, 'PaperUnits','inches','PaperPosition', [0 0 8 8]);

      figure
%      plot(log10(abs(ez_txy(:,4,i))));
      plot(dez08_txy(:,4,i+1));
      hold on
      plot(dez008_txy(:,4,i*10+1),'green');     
      set(gca,'YDir','normal');
%      title(strjoin({'log_{10}[dE_z(x,t=',num2str(frame_index),')]'}));
      title(strjoin({'dE_z(x,t=',num2str(frame_index),')'}));
      legend('\Delta t = 0.08','\Delta t = 0.008','Location','northwest')
      xlim([1 512])
      ylim([-1e-4 1e-4])
      set(gca,'YTick',[-1e-4 -0.5e-4 0 0.5e-4 1e-4]);

%      figure('visible','off');
%      imagesc(log10(abs(dez_txy(:,:,i)))),colorbar;
%      set(gca,'YDir','normal');
%      caxis([-7,0]);
%      title(strjoin({'log_{10}[dE_z(x,y,t=',num2str(frame_index),')]'}));
% the x range 21 to 491 is the range where just the regular box is.
%      xlabel(strjoin({'A_r/A_0 =',num2str(max(squeeze(ez512_txy(21:491,256,i)))/maxv)}))

      saveas(gcf,[num2str(frame_index,'%06i') '.jpg']);
      close gcf;

end
