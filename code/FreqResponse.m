function FreqResponse(Model_obj, Omega, node, ax, linetp )
if nargin < 5
    if nargin == 2
        node = 1;
    elseif nargin < 2
        error('Not enough input arguments.\n Provide: Omega(required), node #s(optional), axes(optional), linetype(optional)',class(Model_obj))
    end
    linetp = '-';
    if nargin < 4
        figure
        ax = axes;
    end
end

for ii = 1:1:length(Omega)
    w = Omega(ii)/60*2*pi;
    Fnew = w^2.*(Model_obj.F);
    Mnew = Model_obj.M;
    Cnew = real(Model_obj.C) + w.*imag(Model_obj.C);
    Knew = real(Model_obj.K) + w.*imag(Model_obj.K);
    
    Z = Knew + 1i*w.*Cnew - w^2.*Mnew;
    X(ii,:) = Z^-1*Fnew;
end
% Hresp(:,1,:) = X';
%     sys = frd(Hresp,Omega/60*2*pi);
hold on
for jj = 1:1:length(node)
    plot(ax, Omega,mean(abs(X(:,node(jj)*4-3:node(jj)*4-2)),2),linetp)
%     index = node(jj)*4-3:node(jj)*4-2;
%     bode(sys(index))
%     subplot(211); hold on
%     plot(Omega,abs(X(:,node(jj)*4-3)),linetp)
%     plot(Omega,abs(X(:,node(jj)*4-2)),linetp)
    ax = gca;
    ax.YScale = 'log';
%     subplot(212); hold on
%     plot(Omega,unwrap(angle(X(:,node(jj)*4-3))),linetp)
%     plot(Omega,unwrap(angle(X(:,node(jj)*4-2))),linetp)
end
hold off
% plot(Omega,max(abs(X)'))
% ax = gca;
% ax.YScale = 'log';