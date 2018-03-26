function FreqResponse(Model_obj, Omega, node, plottype, ax, linetp )
if nargin < 6
    if nargin == 2
        node = 1;
    elseif nargin < 2
        error('Not enough input arguments.\n Provide: Omega(required), node #s(optional), axes(optional), linetype(optional)',class(Model_obj))
    end
    linetp = '-';
    if nargin < 5
        figure
        ax = axes;
    end
    if nargin < 4
        plottype = 'Bode';
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
% Hresp(:,1,:) = X';
%     sys = frd(Hresp,Omega/60*2*pi);
switch plottype
    case 'FreqResponse'
        hold on
        for jj = 1:1:length(node)
            plot(ax, Omega,(abs(X(:,node(jj)*4-2))),linetp)
            ax = gca;
            ax.YScale = 'log';
        end
        hold off
        ax.XLabel.String='Spin Speed \Omega[RPM]';
        ax.YLabel.String='Amplitude[m]';
    case 'Bode'
        hold on
        for jj = 1:1:length(node)
            %     index = node(jj)*4-3:node(jj)*4-2;
            %     bode(sys(index))
            subplot(211); hold on
            plot(Omega,abs(X(:,node(jj)*4-3)),linetp)
            %             plot(Omega,abs(X(:,node(jj)*4-2)),linetp)
            ax = gca;
            ax.YScale = 'log';
            ax.XLabel.String='Spin Speed \Omega[RPM]';
            ax.YLabel.String='Amplitude[m]';
            subplot(212); hold on
            plot(Omega,unwrap(mod(angle(X(:,node(jj)*4-3)),2*pi)),linetp)
            %             plot(Omega,unwrap(angle(X(:,node(jj)*4-2))),linetp)
            ax=gca;
            ax.XLabel.String='Spin Speed \Omega[RPM]';
            ax.YLabel.String='Phase Angle[Rad]';
            ax.YTick=[-3*pi/2:pi/2:0];
            ax.YTickLabel={'-3\pi/2','-\pi','-\pi/2','0'};

        end
        hold off
        

end
% plot(Omega,max(abs(X)'))
% ax = gca;
% ax.YScale = 'log';