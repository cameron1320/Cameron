function Stability( Model_obj, Omega,ax,linetp )
%STABILITY Summary of this function goes here
%   Detailed explanation goes here
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
if nargin < 4
    linetp = '-';
    if nargin < 3
        figure
        ax = axes;
    elseif nargin < 2
        return
    end
    
end

for ii = 1:1:length(Omega)
    w = Omega(ii)/60*2*pi;
    Mnew = Model_obj.M;
    Cnew = real(Model_obj.C) + w.*imag(Model_obj.C);
    Knew = real(Model_obj.K) + w.*imag(Model_obj.K);
    
    Zer = zeros(size(Mnew));
    ey = eye(size(Mnew));
    A = [-Mnew^-1*Cnew, -Mnew^-1*Knew;
        ey,Zer];
    v1=eig(A);
    eiv(:,ii) = sort(v1);
end
while min(abs(imag(eiv(1,:)))) == 0
    eiv(1:2,:) = [];
end
eiv = eiv( ~any( isnan( eiv ) | isinf( eiv ), 2 ),: );
hold on
plot(Omega,max(real(eiv)),'LineWidth',1,'LineStyle',linetp); 
ZeroCross = Omega(zci(max(real(eiv))));
if ~isempty(ZeroCross)
    strcross = ['  Threshold: ' num2str(ZeroCross(1),4) 'RPM'];
    text(ZeroCross(1),0,strcross,'HorizontalAlignment','left','VerticalAlignment','middle','FontWeight','bold','Rotation',90);
    plot(ZeroCross(1),0,'*k'); hold off
end
ax.XAxisLocation = 'origin';
xlabel('Speed (RPM)')
ylabel('Maximum \Re(s)')
title('Stability region')
end

