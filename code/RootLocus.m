function RootLocus(Model_obj,Omega,plotmodes, ax, linetp,color)
%CAMPBELL Summary of this function goes here
%   Detailed explanation goes here
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
if nargin < 6
    color = [0,0,0];
    if nargin <5
    linetp = '.';
    end
    if nargin < 4
        figure
        ax = axes;
    end
    if nargin < 3
        plotmodes = 1:2;
    end
    if nargin < 2
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
%     AA=[Cnew Mnew;Mnew Zer];
%     B=[Knew Zer;Zer -Mnew];
%     v1=eig(B,-1i*AA);
    [~, I] = (sort(abs((v1))));
    eiv(:,ii) = v1((I));
end
% while min(abs(imag(eiv(1,:)))) == 0
%     eiv(1:2,:) = [];
% end
% eiv = eiv( ~any( isnan( eiv ) | isinf( eiv ), 2 ),: );

hold on
for jj = 1:1:length(plotmodes)
   
%     plot(ax,real(eiv(plotmodes(jj)*4-1,:)),(imag(eiv(plotmodes(jj)*4-1,:))),linetp,'Color',[0.8500,0.3250,0.0980])
%     plot(ax,real(eiv(plotmodes(jj)*4-3,:)),(imag(eiv(plotmodes(jj)*4-3,:))),linetp,'Color',[0,0.4470,0.7410])
    plot3(real(eiv(plotmodes(jj)*4-1,:)),abs(imag(eiv(plotmodes(jj)*4,:))),Omega,linetp,'Color',color,'MarkerSize',1);%,'o','Color',[0.8500,0.3250,0.0980],'MarkerSize',3)
    plot3(real(eiv(plotmodes(jj)*4-3,:)),abs(imag(eiv(plotmodes(jj)*4-2,:))),Omega,linetp,'Color',color,'MarkerSize',1);%,'.','Color',[0,0.4470,0.7410],'MarkerSize',4)
%     plot3(real(eiv(plotmodes(jj)*4-0,:)),imag(eiv(plotmodes(jj)*4-0,:)),Omega,'.','Color',[0,0.4470,0.7410])
%     plot3(real(eiv(plotmodes(jj)*4-2,:)),imag(eiv(plotmodes(jj)*4-2,:)),Omega,'.','Color',[0,0.4470,0.7410])
%     r(1,:) = real(eiv(plotmodes(jj)*4-1,:));
%     r(2,:) = real(eiv(plotmodes(jj)*4-3,:));
%     ZeroCross1 = zci(r(1,:));
%     ZeroCross2 = zci(r(2,:));
%     if isempty(ZeroCross1)
%     else
%     strcross1 = [num2str(Omega(ZeroCross1(1)),4) ' (RPM) ' newline];
%     text(ax,0,abs(imag(eiv(plotmodes(jj)*4-1,ZeroCross1(1)))),strcross1,'HorizontalAlignment','right','VerticalAlignment','bottom','FontWeight','bold');
%     plot(ax,0,abs(imag(eiv(plotmodes(jj)*4-1,ZeroCross1(1)))),'*k');
%     end
%     if isempty(ZeroCross2)
%     else
%     strcross2 = [ num2str(Omega(ZeroCross2(1)),4) ' (RPM) ' newline];
%     text(ax,0,abs(imag(eiv(plotmodes(jj)*4-3,ZeroCross2(1)))),strcross2,'HorizontalAlignment','right','VerticalAlignment','bottom','FontWeight','bold');
%     plot(ax,0,abs(imag(eiv(plotmodes(jj)*4-3,ZeroCross2(1)))),'*k');
%     end
end
hold off
ax.YAxisLocation = 'origin';
axis([-inf,inf,0,inf])
xlabel('Real(s)')
ylabel('Imag(s)')
end
