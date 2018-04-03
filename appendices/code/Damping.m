function Damping(Model_obj,Omega,plotmodes,ax,linetp,ch)
%CAMPBELL Summary of this function goes here
%   Detailed explanation goes here
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);

if nargin < 6
    if nargin == 2
        plotmodes = 1:2;
    elseif nargin < 2
        return
    end
    ch = 'dec';
    if nargin < 5
        linetp = '-';
    end
    if nargin < 4
        figure
        ax = axes;
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
        [~, I] = (sort(abs((v1))));
    eiv(:,ii) = v1((I));
end

while min(abs(imag(eiv(1,:)))) == 0
    eiv(1:2,:) = [];
end
eiv = eiv( ~any( isnan( eiv ) | isinf( eiv ), 2 ),: );

hold on
for jj = 1:1:length(plotmodes)
    idf = plotmodes(jj)*4-1; %forward mode index location
    idb = idf - 2; %backward mode index location
    zeta = -real(eiv([idb,idf],:))./abs(eiv([idb,idf],:));
%     delta = 2*pi.*zeta./sqrt(1-zeta.^2); %Logarithmic Decriment
%     delta = -2*pi*real(eiv([idb,idf],:))./imag(eiv([idb,idf],:));
    plot(Omega,zeta(1,:),linetp,'Color',[0,0.4470,0.7410],'LineWidth',1);
    plot(Omega,zeta(2,:),linetp,'Color',[0.8500,0.3250,0.0980],'LineWidth',1);
    ZeroCross1 = Omega(zci(zeta(1,:)));
    ZeroCross2 = Omega(zci(zeta(2,:)));
    if isempty(ZeroCross1)
    else
    strcross1 = [num2str(ZeroCross1(1),4) newline];
    text(ZeroCross1(1),0,strcross1,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold');
    plot(ZeroCross1(1),0,'*k');
    end
    if isempty(ZeroCross2)
    else
    strcross2 = [ num2str(ZeroCross2(1),4) newline];
    text(ZeroCross2(1),0,strcross2,'HorizontalAlignment','left','VerticalAlignment','bottom','FontWeight','bold');
    plot(ZeroCross2(1),0,'*k');
    end

end
hold off
ax.XAxisLocation = 'origin';
xlabel('Speed (RPM)')
ylabel('\zeta')
title('Damping Characteristics')

end
