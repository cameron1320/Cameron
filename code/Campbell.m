function Campbell(Model_obj,Omega,plotmodes,ax,linetp )
%CAMPBELL Summary of this function goes here
%   Detailed explanation goes here
if nargin < 5
    if nargin == 2
        plotmodes = 1:2;
    elseif nargin < 2
        return
    end
    linetp = '.';
    if nargin < 4
        ax = [];
    end
end
for ii = 1:1:length(Omega)
    w = Omega(ii)/60*2*pi;
    Mnew = Model_obj.M;
    Cnew = real(Model_obj.C) + w.*imag(Model_obj.C);
    Knew = real(Model_obj.K) + w.*imag(Model_obj.K);
    
    Zer = zeros(size(Mnew));
    AA=[Cnew Mnew;Mnew Zer];
    B=[Knew Zer;Zer -Mnew];
    v1=eig(B,-AA);
    [~, I] = (sort(abs(v1)));
    eiv(:,ii) = v1(I);
end
if isempty(ax)
    figure
    ax = axes;
end
% while min(abs(imag(eiv(1,:)))) == 0
%     eiv(1:2,:) = [];
% end
% eiv = eiv( ~any( isnan( eiv ) | isinf( eiv ), 2 ),: );

hold on
for jj = 1:1:length(plotmodes)
    plot(ax,Omega,abs(imag(eiv(plotmodes(jj)*4-3,:)))/2/pi*60,linetp,'Color',[0,0.4470,0.7410],'LineWidth',1)
    plot(ax,Omega,abs(imag(eiv(plotmodes(jj)*4-1,:)))/2/pi*60,linetp, 'Color',[0.8500,0.3250,0.0980],'LineWidth',1)
%         plot(ax,Omega,imag(eiv(plotmodes(jj)*4-2,:))/2/pi*60,linetp,'Color',[0,0.4470,0.7410],'LineWidth',1)
%     plot(ax,Omega,imag(eiv(plotmodes(jj)*4,:))/2/pi*60,linetp, 'Color',[0.8500,0.3250,0.0980],'LineWidth',1)

%     plot(ax,Omega,abs(eiv(plotmodes(jj),:))*60/2/pi,'Color',[0,0.4470,0.7410],'LineWidth',1,'LineStyle',linetp)
    
end
plot(Omega,Omega,'--k','LineWidth',2)
hold off
xlabel('Spin speed (RPM)')
ylabel('Whirl speed (RPM)')
title('Campbell Diagram')


end

