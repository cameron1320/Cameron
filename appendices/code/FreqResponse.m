function FreqResponse(Model_obj, Omega, node, plottype, fi, linetp )
if nargin < 6
    if nargin == 2
        node = 1;
    elseif nargin < 2
        error('Not enough input arguments.\n Provide: Omega(required), node #s(optional), axes(optional), linetype(optional)',class(Model_obj))
    end
    linetp = '-';
    if nargin < 5
        fi = figure;
        ax = axes(fi);
        if nargin < 4
            plottype = 'Bode';
        end
        if strcmp(plottype,'Bode') == 1
            ax(1)=subplot(2,1,1,ax);
            ax(2)=subplot(2,1,2);
        end
    end

end

for ii = 1:1:length(Omega)
    w = Omega(ii)/60*2*pi;
    Fnew = w^2.*(Model_obj.F);
    Mnew = Model_obj.M;
    Cnew = real(Model_obj.C) + w.*imag(Model_obj.C);
    Knew = real(Model_obj.K) + w.*imag(Model_obj.K);
    n=length(Mnew);
    I=kron(eye(n/2),diag([1,1i]));
    Z = I*(Knew + 1i*w.*Cnew - w^2.*Mnew);
    X(ii,:) = Z^-1*Fnew;
end
switch plottype
    case 'FreqResponse'
        hold on
        for jj = 1:1:length(node)
            ax=fi.Children;
            plot(ax, Omega,(abs(X(:,node(jj)*4-3))),linetp)
            plot(ax, Omega,(abs(X(:,node(jj)*4-2))),linetp)
            ax.YScale = 'log';
        end
        hold off
        ax.XLabel.String='Spin Speed \Omega[RPM]';
        ax.YLabel.String='Amplitude[m]';
    case 'Bode'
        hold on
        for jj = 1:1:length(node)
            subplot(211); hold on
            ax(1)=fi.Children(1);
%             plot(Omega,abs(X(:,node(jj)*4-3)),linetp)
%             plot(Omega,abs(X(:,node(jj)*4-2)),linetp)
            plot(ax(1),Omega,abs(39370.1*X(:,node(jj)*4-3)),linetp)
            plot(ax(1),Omega,abs(39370.1*X(:,node(jj)*4-2)),linetp)
            ax(1).YScale = 'log';
            ax(1).XLabel.String='Spin Speed \Omega[RPM]';
%             ax(1).YLabel.String='Amplitude[m]';
            ax(1).YLabel.String='Amplitude[mils]';
            subplot(212); hold on
            ax(2)=fi.Children(2);
%             plot(Omega,unwrap(mod(angle(X(:,node(jj)*4-3)),2*pi)),linetp)
%             plot(Omega,unwrap(mod(angle(X(:,node(jj)*4-2)),2*pi)),linetp)
            plot(ax(2),Omega,180/pi*unwrap(mod(angle(X(:,node(jj)*4-3)),2*pi)),linetp)
            plot(ax(2),Omega,180/pi*unwrap(mod(angle(X(:,node(jj)*4-2)),2*pi)),linetp)
            ax(2).XLabel.String='Spin Speed \Omega[RPM]';
%             ax.YLabel.String='Phase Angle[Rad]';
            ax(2).YLabel.String='Phase Angle[deg.]';
%             ax.YTick=[-3*pi/2:pi/2:0];
%             ax.YTickLabel={'-3\pi/2','-\pi','-\pi/2','0'};
            ax(2).YTick=[-360:90:0];
            ax(2).YTickLabel={'0','-270','-180','-90','0'};
        end
        hold off
        

end
% plot(Omega,max(abs(X)'))
% ax = gca;
% ax.YScale = 'log';