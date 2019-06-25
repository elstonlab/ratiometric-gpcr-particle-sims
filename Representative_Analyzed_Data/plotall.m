function plotall()
% Loads processeed data and plots all.
% The curves correspond to results shown in Figure 5G, 6B, and 8E.
n10000_r_unif=csvread('NR10000_uniform_ratiometric.csv');
n10000_nr_unif=csvread('NR10000_uniform_nonratiometric.csv');
n10000_r_align=csvread('NR10000_align_ratiometric.csv');
n10000_nr_align=csvread('NR10000_align_nonratiometric.csv');
n10000_r_orthogonal=csvread('NR10000_orthogonal_ratiometric.csv');
n10000_nr_orthogonal=csvread('NR10000_orthogonal_nonratiometric.csv');
n10000_r_diagonal=csvread('NR10000_diagonal_ratiometric.csv');
n10000_nr_diagonal=csvread('NR10000_diagonal_nonratiometric.csv');
n10000_r_oppose=csvread('NR10000_oppose_ratiometric.csv');
n10000_nr_oppose=csvread('NR10000_oppose_nonratiometric.csv');

fh=figure('position',[234 1024 1515 251]);
addsubplot(fh,1,n10000_nr_unif,n10000_r_unif,'uniform');
addsubplot(fh,2,n10000_nr_align,n10000_r_align,'align');
addsubplot(fh,3,n10000_nr_orthogonal,n10000_r_orthogonal,'orthogonal');
addsubplot(fh,4,n10000_nr_diagonal,n10000_r_diagonal,'diagonal');
addsubplot(fh,5,n10000_nr_oppose,n10000_r_oppose,'oppose');

savefig('angle_over_time');
print(gcf,'-dtiff','angle_over_time.tif','-r300');
end

function addsubplot(fighandle,subplotidx,nr_data,r_data,subplottitle)
figure(fighandle);
subplot(1,5,subplotidx); hold on;
if strcmpi(subplottitle,'align') % Flip order to show both sets of curves.
    a=plot(r_data(:,1),r_data(:,2:end),'color',[0.8500    0.3250    0.0980]);    
    b=plot(nr_data(:,1),nr_data(:,2:end),'color',[0 .447 .741]);
else
    b=plot(nr_data(:,1),nr_data(:,2:end),'color',[0 .447 .741]);
    a=plot(r_data(:,1),r_data(:,2:end),'color',[0.8500    0.3250    0.0980]);
end
xlabel('time (s)');
ylabel('angle (degrees)')
legend([a(1),b(1)],'Ratiometric','Nonratiometric','location','northeast');
set(gca,'xlim',[0 600],'xtick',0:100:600,'ylim',[-5 185],'ytick',0:30:180);
set(gca,'fontsize',14)
title(subplottitle)
end