function varargout = baloa( ...
    x, xName, y, yName, a, ...
    doPlotMD, axMD, ...
    doPlotC, axC, ...
    doPlotBasicStats, doPlotExtStats, doPlotLS)
%% preparation
% multiFig = false;
if doPlotMD || doPlotC
    % if doMDPlot, f(1) = axMD.Parent; end
    % if doCPlot, f(2) = axC.Parent; end
    ax = [axMD;axC];
    f = get(ax,'Parent');
    if iscell(f), f = vertcat(f{:}); end
    f = unique(f);
    % f can be the handle to one or more figures
    % if numel(f)>1, multiFig = true; end
else
    f = [];
end

% keep only values that can be used in calculations
lok = isfinite(x) & isnumeric(x) & isfinite(y) & isnumeric(y);
n = nnz(lok);
xok = x(lok);
yok = y(lok);

%% calculation
% mean statistics
muXY = mean([xok,yok],2);

% difference statistics
d = xok-yok; % difference
muD = mean(d); % mean of difference

sD = std(d); % s_d in article
% varD = var(d); % s_d^2 in article

varMuD = sD^2/n; % variance of muD (SE^2) (p. 141)
seMuD = sqrt(varMuD); % standard error of muD (p. 142)

% varSD = sD^2/2/(n-1); % approximated variance of sD (article p. 141)
% Instead, use exact formula for unbiased estimated variance
% Source: http://stats.stackexchange.com/a/28567/80486
% sSD = sD * gamma((n-1)/2) / gamma(n/2) * ...
%     sqrt( (n-1)/2 - ( gamma(n/2) / gamma((n-1)/2) )^2 );
gammafrac = gamma((n-1)/2) / gamma(n/2);
% if ~isfinite(gammafrac) % true for large n
% % approximate using gamma(a+b)/gamma(a) ~ a^b
% % Source: https://en.wikipedia.org/w/index.php?title=Gamma_function&oldid=737220343#General
% % compare:
% % figure
% % n = 1:500;
% % g1 = gamma((n-1)/2)./gamma(n/2);
% % g2 = ((n-1)/2).^(-1/2);
% % g3 = sqrt(2./(n-1));
% % plot(n,g1,n,g2,n,g3,n,g1-g2,n,g1-g3,n,g2-g3)
% % legend g1 g2 g3 g1-g2 g1-g3 g2-g3
% gammafrac = sqrt(2/(n-1)); % same as: gammafrac = ((n-1)/2).^(-1/2);
% end
sSD = sD * gammafrac * sqrt( (n-1)/2 - gammafrac^-2 );
varSD = sSD^2; % unbiased estimate of variance of s_d

% limits of agreement statistics
p = 1-a/2;
z = Ninv(p); % inverse normal distribution at p = 1-alpha/2
loa = muD + z*sD*[-1 1]; % limits of agreement (LOA)

% confidence intervals (CI) for muD and loa
t = Tinv(p,n-1);
varLoa = varMuD + z^2*varSD; % article: bottom of p. 141
seLoa = sqrt(varLoa); % standard error of the LOA
eLoa = t*seLoa; % LOA error
eMuD = t*seMuD; % muD error
muDCI = muD + eMuD*[-1 1];
loaCI = [loa;loa] + eLoa*[-1 -1;1 1];
% loaCI is a 2x2 matrix. Every column contains the CI for a LOA: the first
% column corresponds to loa(1), the second to loa(2). The first and second
% row correspond to the lower and upper CI bound respectively.
% LOA = [loaCI(1,:);loa;loaCI(2,:)]; % optional 3x2 matrix form

% mean difference correlation statistics
[rSMuD,pRSMuD] = corr(muXY,d,'type','Spearman'); %TODO make independent of stats toolbox?

% linear regression and correlation %TODO add polyfit for muXY and d
[polyXY,statsPXY] = polyfit(xok,yok,1);
[rhoXY,pRhoXY] = corrcoef(xok,yok);
rhoXY = rhoXY(1,2);
pRhoXY = pRhoXY(1,2);

% linear regression statistics
resPXY = yok - polyval(polyXY,xok); % residual
ssePXY = sum(resPXY.^2); % or rss/ssr
msePXY = ssePXY/statsPXY.df;
% compare with the same calculation (CF toolbox required):
% [~,gof] = fit(xok,yok,'Poly1');
% gof.sse-sse % equal within double precision
% gof.rmse-sqrt(mse) % idem
% gof.dfe-s.df % equal

%% graphics
% correlation plot
if doPlotC
    % preparation
    axes(axC)
    legEntries = gobjects(0);
    
    % plot y against x
    sC = scatter(xok,yok);
    sC.UserData = dcStruct([],'M1','M2',[],@dcXY);
    axis tight
    axis equal
    
    if doPlotBasicStats
        % add correlation to legend
        if pRhoXY<1e-4
            % as recommended by BMJ 1996;312:572
            % http://dx.doi.org/10.1136/bmj.312.7030.572
            strPRhoXY = sprintf('\\itp\\rm < 0.0001');
        else
            strPRhoXY = sprintf('\\itp\\rm = %.2f',pRhoXY);
        end
        sC.DisplayName = sprintf('\\rho = %.2f (%s)',rhoXY,strPRhoXY);
        legEntries(end+1) = sC;
    end
    
    % plot least-squares line
    if doPlotLS
        XYLine = refline;
        XYLine.Color = 'r';
        legEntries(end+1) = XYLine;
        XYLine.DisplayName = 'least-squares';
        XYLine.UserData = dcStruct( ...
            ['M2 = ' num2str(polyXY(1)) 'M1 + ' num2str(polyXY(2))], ...
            [], [], ...
            ['MSE = ' num2str(msePXY)], ...
            @dcXY);
    end
    
    % y = x reference line
    eqLine = refline(1,0); % 45 degree, because of axis equal
    eqLine.Color = [.75 .75 .75];
    eqLine.LineStyle = '--';
    eqLine.DisplayName = 'line of equality';
    eqLine.UserData = dcStruct([],[],[],'M2 = M1',@dcXY);
    legEntries(end+1) = eqLine;
    
    % axes labels
    xlabel('M_1')
    ylabel('M_2')
    title(sprintf(['Scatter plot of (%u observations):\n' ...
        ' \\rm\\itM_1\\rm: %s\n \\itM_2\\rm: %s'],n,xName,yName))
    
    % legend
    legend(legEntries,'Location','SouthEast')
    
    % reorder plot children
    axC.Children = axC.Children([end 1:end-1]);
end

% mean-difference plot
if doPlotMD
    % preparation
    axes(axMD)
    legEntries = gobjects(0);
    
    % mean-difference plot
    sMD = scatter(muXY,d);
    sMD.UserData = dcStruct([],'?','d',[],@dcXY); %TODO check ?
    axis equal
    xl = xlim;
    
    % plot statistics
    if doPlotBasicStats
        % prepare axes
        hold(axMD,'on')
        padding = range(ylim)/20; % distance to enhance visibility
        yl = ylim;
        yl = [ ...
            min(yl(1), loaCI(1,1)) - padding, ...
            max(yl(2), loaCI(2,2)) + padding ...
            ];
        ylim(yl)
        
        % add correlation to legend
        if pRSMuD<1e-4
            % as recommended by BMJ 1996;312:572
            % http://dx.doi.org/10.1136/bmj.312.7030.572
            strPRSMuD = sprintf('\\itp\\rm < 0.0001');
        else
            strPRSMuD = sprintf('\\itp\\rm = %.2f',pRSMuD);
        end
        sMD.DisplayName = sprintf( ...
            '\\itr_s\\rm = %.2f (%s)',rSMuD,strPRSMuD);
        legEntries(end+1) = sMD;
        
        % lower LOA line
        lLine = refline(0,loa(1));
        lLine.Color = 'k';
        lLine.UserData = dcStruct([],[],'lower LOA', ...
            [num2str(100*(1-a)) '% CI = [' num2str(loaCI(1,1)) ', ' num2str(loaCI(2,1)) ']'], ...
            @dcXY);
        text(xl(2)-padding/2,loa(1)+padding/2,sprintf( ...
            '$\\overline{d}-%.2fs_d$',z ...
            ), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        
        % upper limit of agreement line
        uLine = refline(0,loa(2));
        uLine.Color = 'k';
        uLine.UserData = dcStruct([],[],'upper LOA', ...
            [num2str(100*(1-a)) '% CI = [' ...
            num2str(loaCI(1,2)) ', ' num2str(loaCI(2,2)) ']'], ...
            @dcXY);
        text(xl(2)-padding/2,loa(2)+padding/2, ...
            sprintf( ...
            '$\\overline{d}+%.2fs_d$',z ...
            ), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        
        % mean difference line
        muDLine = refline(0,muD); % mean difference
        muDLine.Color = 'k';
        text(xl(2)-padding/2,muD+padding/2,'$\overline{d}$', ...
            'Interpreter','latex', ...
            'HorizontalAlignment','right')
        muDLine.UserData = dcStruct([],[],'mean difference', ...
            [num2str(100*(1-a)) '% CI = [' ...
            num2str(muDCI(1)) ', ' num2str(muDCI(2)) ']'], ...
            @dcXY);
        
        % plot additional statistics
        if doPlotExtStats
            % lower LOA errorbar
            lE = errorbar(xl(2),loa(1),eLoa,'r');
            lE.UserData = lLine.UserData;
            
            % mean difference errorbar
            muDE = errorbar(xl(2),muD,eMuD,'r');
            muDE.UserData = muDLine.UserData;
            
            % upper LOA errorbar
            uE = errorbar(xl(2),loa(2),eLoa,'r');
            uE.UserData = uLine.UserData;
        end
    end
    
    % plot least-squares line
    if doPlotLS
        % least-squares line
        lsMuDLine = refline;
        lsMuDLine.Color = 'r';
        lsMuDLine.DisplayName = ...
            'least-squares';
        lsMuDLine.UserData = dcStruct( ...
            'least-squares line of mean and difference', ...
            [], [], [], ... %TODO add parameters from polyfit like in correlation plot
            @dcXY);
        legEntries(end+1) = lsMuDLine;
    end
    
    % zero difference line
    line0 = refline(0,0);
    line0.Color = [.75 .75 .75];
    line0.LineStyle = '--';
    line0.DisplayName = 'line of equality';
    line0.UserData = dcStruct(line0.DisplayName,[],[],[],@dcXY);
    legEntries(end+1) = line0;
    
    % axes labels
    xlabel(sprintf('mean \\it = (M_1+M_2)/2')) %TODO check ?
    ylabel(sprintf('difference \\itd = M_1-M_2'))
    title(sprintf(['Mean-difference plot of (%u observation pairs):\n' ...
        ' \\rm\\itM_1\\rm: %s\n \\itM_2\\rm: %s'],n,xName,yName))
    
    % legend
    legend(legEntries,'Location','SouthWest')
    
    % reorder plot children
    axMD.Children = axMD.Children([end 1:end-1]);
end

%% set data cursor update function for figure(s)
for f = f(:).'
    dc = datacursormode(f);
    dc.UpdateFcn = @dcUpdateFcn;
    dc.SnapToDataVertex = 'off';
    dc.Enable = 'on';
end

%% output
out.muD = muD;
out.muDCI = muDCI;
out.loa = loa;
out.loaCI = loaCI;
out.sD = sD;
out.rSMuD = rSMuD;
out.pRSMuD = pRSMuD;
varargout = {{out}};
end