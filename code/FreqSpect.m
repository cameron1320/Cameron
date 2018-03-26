function FreqSpect(Model_obj, Omega, range, node, ax, linetp )
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
w = range./60.*2.*pi;
for ii = 1:1:length(Omega)
    for jj = 1:1:length(w)
        O = Omega(ii)/60*2*pi;
        Fnew = O^2.*(Model_obj.F);
        Mnew = Model_obj.M;
        Cnew = real(Model_obj.C) + 1i*O.*imag(Model_obj.C);
        Knew = real(Model_obj.K) + 1i*O.*imag(Model_obj.K);
        
        Z = Knew + 1i*w(jj).*Cnew - w(jj)^2.*Mnew;
        X(ii,jj,:) = Z^-1*Fnew;
    end
end
% Hresp(:,1,:) = X';
%     sys = frd(Hresp,Omega/60*2*pi);
hold on
% Y=squeeze(X(:,100,:));
% [W,Om] = meshgrid(range,Omega);
for jj = 1:1:length(node)
P=squeeze(mean(abs(X(:,:,node(jj)*4-3:node(jj)*4-2)),3));
whos;
%     waterfall(ax, Omega,w,X(:,:,node(jj)*4),linetp)
    waterfall(range,Omega,P)
%     index = node(jj)*4-3:node(jj)*4-2;
%     bode(sys(index))
%     subplot(211); hold on
%     plot(Omega,abs(X(:,node(jj)*4-3)),linetp)
%     plot(Omega,abs(X(:,node(jj)*4-2)),linetp)
    ax = gca;
    ax.YScale = 'linear';
%     subplot(212); hold on
%     plot(Omega,unwrap(angle(X(:,node(jj)*4-3))),linetp)
%     plot(Omega,unwrap(angle(X(:,node(jj)*4-2))),linetp)
end
hold off
% plot(Omega,max(abs(X)'))
ax = gca;
ax.ZScale = 'log';