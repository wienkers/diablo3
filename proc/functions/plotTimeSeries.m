function f = plotTimeSeries(time,gzf,dat,name,normalise)

if nargin < 5
    normalise = false;
end
if nargin < 4
    name = '';
end

if size(dat,1) ~= length(gzf)
    dat = dat.';
end


rdbu = flipud(cbrewer('div', 'RdBu', 128));

zmask = gzf > 0.1*max(gzf) & gzf < 0.9*max(gzf);
databs = max(abs(dat(zmask,:)),[],1);

tmask = time > 1;

databs = smooth(databs,100).';

% if length(gzf) > 1024
%     skip = fix(length(gzf)/1024);
%     dat = dat(1:skip:end,:);
%     gzf = gzf(1:skip:end);
%     databs = databs(1:skip:end,:);
% end

if normalise

    dat = dat./databs;
    n_panels = 5;%6;

    f = figure();
    ax(1) = subplot(n_panels,1,1:4);
    pcolor(time,gzf,dat); shading interp;
    colormap(rdbu);
    caxis([-1,1]);
    title(name)
    ylabel('$z$')

    ax(2) = subplot(n_panels,1,5);
    semilogy(time,databs);
    ylabel('Magnitude')
    if max(databs(:,tmask)) ~= 0
        ylim([min(databs(:,tmask)),1.5*max(databs(:,tmask))]);
    end


%     ax(3) = subplot(n_panels,1,n_panels);
%     semilogy(time,tke_int);
%     ylabel('$TKE$')


    set(ax(1),'XTickLabel','');
    set(ax(1:2),'layer','top');
    xlabel('$t \; [f^{-1}]$')
    linkaxes(ax,'x');
    pause(0.2);
    xlim([min(time),max(time)]);
    pause(0.2);
    axes(ax(1));

else

    f = figure();
    pcolor(time,gzf,dat); shading interp;
    colourbar(name,rdbu, max(databs)*[-1,1]);
    ylabel('$z$')
    xlabel('$t$')

end



end
