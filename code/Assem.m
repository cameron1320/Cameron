function Assem(obj,input_file)
run(input_file)

nodes = size(elems,1)+1;
npos = zeros(nodes,1);
for ii = 1:nodes-1
    npos(ii+1) = npos(ii) + elems(ii,2);
end
obj.npos = npos;

K = zeros(nodes*DOF,nodes*DOF);
C = zeros(nodes*DOF,nodes*DOF);
M = zeros(nodes*DOF,nodes*DOF);
F = zeros(nodes*DOF,1);
for ii = 1:1:size(elems,1)
    E = mate(elems(ii,1),2);
    l = elems(ii,2);
    do = elems(ii,3);
    if size(elems(ii,:))<4
        di=0;
    else
        do=elems(ii,4);
    end
    I = obj.AInert(do,di);
    A = obj.Area(do,di);
    rho = mate(elems(ii,1),3);
    nuv = mate(elems(ii,1),4);
    poi = mate(elems(ii,1),5);
    Id = obj.DiaInert(do,l,rho,di);
    Ip = obj.PolInert(do,l,rho,di);
    [kbe, cbe, mbe] = obj.TBeam(E,I,l,A,rho,Id,Ip,nuv,poi,DOF);
    K = obj.Add(K,kbe,(ii-1)*DOF+1:(ii+1)*DOF);
    C = obj.Add(C,cbe,(ii-1)*DOF+1:(ii+1)*DOF);
    M = obj.Add(M,mbe,(ii-1)*DOF+1:(ii+1)*DOF);
end

if exist('disks', 'var')
    for ii = 1:1:size(disks,1)
        do = disks(ii,3);
        l = disks(ii,2);
        rho = mate(disks(ii,1),3);
        a = disks(ii,4);
        Ki = disks(ii,5);
        if size(disks(ii,:))<7
            di=0;
        else
            di=disks(ii,7);
        end
        md = obj.Mass(do,l,rho,di);
        Ip = obj.PolInert(do,l,rho,di);
        Id = obj.DiaInert(do,l,rho,di);
        [mi,gi,fi] = obj.Disk(md,Id,Ip,a,Ki,DOF);
        M = obj.Add(M,mi,(disks(ii,6)*DOF-DOF+1:disks(ii,6)*DOF));
        C = obj.Add(C,gi,(disks(ii,6)*DOF-DOF+1:disks(ii,6)*DOF));
        F = obj.AddVect(F,fi,(disks(ii,6)*DOF-DOF+1:disks(ii,6)*DOF));
    end
end

for ii = 1:1:size(bears,1)
    kx = beartypes(bears(ii,1),2);
    ky = beartypes(bears(ii,1),3);
    c = beartypes(bears(ii,1),4);
    [ki,ci] = obj.Bear(kx,ky,c,DOF);
    K = obj.Add(K,ki,(bears(ii,2)*DOF-DOF+1:bears(ii,2)*DOF));
    C = obj.Add(C,ci,(bears(ii,2)*DOF-DOF+1:bears(ii,2)*DOF));
    
end
if exist('mags', 'var')
    for ii = 1:1:size(mags,1)
        olddir = pwd;
        cd(char(magfile(mags(ii,1),2)))
        Magfn = str2func(char(magfile(mags(ii,1),3)));
        [ki, ci] = Magfn();
        cd(olddir)
        K = obj.Add(K,ki,(mags(ii,2)*DOF-DOF+1:mags(ii,2)*DOF));
        C = obj.Add(C,ci,(mags(ii,2)*DOF-DOF+1:mags(ii,2)*DOF));
        
    end
end
obj.M = M;
obj.C = C;
obj.K = K;
obj.F = F;