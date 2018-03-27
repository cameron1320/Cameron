function Shape(Model_obj,Omega,plotmodes,ax,linetp )
%CAMPBELL Summary of this function goes here
%   Detailed explanation goes here
if nargin < 5
    if nargin < 2
        return
    end
    linetp = '-';
    if nargin < 4
        figure
        ax = axes;
    end
end

syms x;
l=1;
z=x/l;
alph=.25;
Tt1 = (1/(1+alph))*(2*z^3 - 3*z^2 - alph*z + 1 + alph);
Tt2 = (1/(1+alph))*(-2*z^3+3*z^2+alph*z);
Tr1 = (l/(1+alph))*(z^3 - (2+1/2*alph)*z^2 + (1+1/2*alph)*z);
Tr2 = (l/(1+alph))*(z^3 - (1-1/2*alph)*z^2 - 1/2*alph*z);
cpsi = [Tt1  0    0  Tr1  Tt2  0    0  Tr2;
       0  Tt1  -Tr1  0    0  Tt2  -Tr2  0];
% psi = [1 - 3*(x/l)^2 + 2*(x/l)^3, x - 2*l*(x/l)^2 + l*(x/l)^3, 3*(x/l)^2 - 2*(x/l)^3, -l*(x/l)^2 + l*(x/l)^3];
% cpsi = [psi(1),      0,       0, psi(2), psi(3),      0,       0, psi(4);
%     0, psi(1), -psi(2),      0,      0, psi(3), -psi(4),      0];
ipts = 4;
d = (0:ipts-1)/ipts;
for ii = 1:1:length(d)
    Nx(ii,:) = subs(cpsi(1,:),x,d(ii));
    Ny(ii,:) = subs(cpsi(2,:),x,d(ii));
end
Nx = double(Nx);
Ny = double(Ny);
w = Omega/60*2*pi;
Mnew = Model_obj.M;
Cnew = real(Model_obj.C) + w.*imag(Model_obj.C);
Knew = real(Model_obj.K) + w.*imag(Model_obj.K);
nM = length(Mnew);

% Zer = zeros(nM);
% AA=[Cnew Mnew;Mnew Zer];
% B=[Knew Zer;Zer -Mnew];
    Zer = zeros(size(Mnew));
    ey = eye(size(Mnew));
    A = [-Mnew^-1*Cnew, -Mnew^-1*Knew;
        ey,Zer];

[V,D] = eig(A,'vector');
% ind=find( (imag(D) > 0) & (abs(D)>1e-3) );	% take positive frequencies only
% l1=D(ind);
% v1=V(nM+1:end,ind);  	
l1=D;
v1=V(nM+1:end,:);
% v1=V(1:nM,:);
[~, I] = sort(abs((l1)));
l2 = D(I);
v2 = V(:,I);
v = v2(:,plotmodes);
v=v/norm(v);
for ii = 1:1:length(Model_obj.M)/4-1
    Li = Model_obj.npos(ii+1)-Model_obj.npos(ii);
    Nxnew = Nx*diag([1 1 1 Li 1 1 1 Li]);
    Nynew = Ny*diag([1 1 Li 1 1 1 Li 1]);
    pos(ipts*(ii-1)+1:ipts*(ii-1)+ipts+1) = linspace(Model_obj.npos(ii),Model_obj.npos(ii+1),ipts+1);
    pos(end)=[];
    shapex(ipts*(ii-1)+1:ipts*(ii-1)+ipts) = Nxnew*v(4*(ii-1)+1:4*(ii-1)+8);
    shapey(ipts*(ii-1)+1:ipts*(ii-1)+ipts) = Nynew*v(4*(ii-1)+1:4*(ii-1)+8);
end
pos = [pos, Model_obj.npos(end)];
shapex = [shapex, v(nM-3)];
shapey = [shapey, v(nM-2)];
plot(ax,[pos;pos]',real([shapex; shapey])')
hold on
for ii = 1:length(Model_obj.npos)
    plot([Model_obj.npos(ii);Model_obj.npos(ii)]',real([shapex(ipts*(ii-1)+1); shapey(ipts*(ii-1)+1)])','.k','MarkerSize',8)
end
hold off
figure
cpts = 20;
theta = linspace(0,2*pi,cpts);
for ii =1:length(shapex)
    xx(:,ii)=real(shapex(ii))*cos(theta) - imag(shapex(ii))*sin(theta);
    y(:,ii)=real(shapey(ii))*cos(theta) - imag(shapey(ii))*sin(theta);
end
N = repmat(pos,cpts,1);
hold on
plot3(real(shapex),pos,real(shapey),'-k')

plot3(xx,N,y,'Color',[0,0.4470,0.7410])

for ii = 1:length(Model_obj.npos)
    plot3(real(shapex(ipts*(ii-1)+1)),Model_obj.npos(ii),real(shapey(ipts*(ii-1)+1)),'.k','MarkerSize',8)
end

hold off
mx = max([abs(shapex),abs(shapey)]);
axis([-mx, mx, -inf, inf,-mx, mx ]);
disp(['Damping: ' num2str(real(l2(plotmodes))) '. Frequency: ' num2str(imag(l2(plotmodes))/2/pi*60) '(RPM).'])
ax=gca;
ax.XDir='reverse';
ax.Title.String={(['\Re(s): ' num2str(real(l2(plotmodes)),2) ', \Im(s): ' num2str(imag(l2(plotmodes))/2/pi*60,4) '[RPM]']);(['\Omega: ' num2str(Omega,4) '[RPM]'])};
% 
% [V, D] = eig(Knew, Mnew,'vector');
% D = sqrt(D)*60/2/pi;
% [D, ind] = sort(abs(D));
% D(plotmodes)
% V = V(:,ind);
% v2 = (V(:,plotmodes));
% v3 = reshape(v2,4,[]);
% v4 = bsxfun(@times, v3, 1./sqrt(sum(v3.^2, 2)));
% v2 = reshape(v4,1,[])';
% for ii = 1:1:length(Model_obj.M)/4-1
%     shape(8*(ii-1)+1:8*(ii-1)+8) = N*v2(4*(ii-1)+1:4*(ii-1)+8);
% end
% shape = [shape, v2(end-3)];
% plot(real(shape))


% if isempty(ax)
%     figure
%     ax = axes;
% end
% 
% hold on
% for jj = 1:1:length(plotmodes)
%     pv(1,:) = real(eivect(speed,[4.*(1:14)-3],plotmodes(jj)));
%     
%     plot(pv)
% end
% hold off


end
