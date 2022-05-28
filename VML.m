clear ; close all; clc;

%% Parametro

% Parametros Ala Principal
prompt = {'What is the chord? [m]','What is the Aspect Ratio?', 'What is the sweep? ', 'What is the taper ratio? ', 'What is the dihedral? ', 'Insert AOA', 'What is the sideslip angle? ', 'What is the mach? ', 'What is the height? ','Number control points', 'Number of panels around span'   };
dlgtitle = 'Wing Parametro';
dims = [1 35];
definput = {'5','5','50','0.3','0','2.1,3','0','0.3','1000','10','20'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
Ala.Parametro.Chord  = str2num(answer{1});
Ala.Parametro.AR = str2num(answer{2});
Ala.Parametro.Sweep = str2num(answer{3});
Ala.Parametro.Sweepr = deg2rad(Ala.Parametro.Sweep);
Ala.Parametro.Sw = Ala.Parametro.AR*Ala.Parametro.Chord^2;
Ala.Parametro.TR = str2num(answer{4});
Ala.Parametro.b = (Ala.Parametro.AR*Ala.Parametro.Sw)^(0.5);
Ala.Parametro.dihedral = str2num(answer{5});
Ala.Parametro.dihedralr = deg2rad(Ala.Parametro.dihedral);


Ala.Parametro.xoffset   = 0;            
Ala.Parametro.yoffset   = 0;            
Ala.Parametro.zoffset   = 0;            

Ala.Parametro.prof_R  = '0000';
Ala.Parametro.prof_T  = '0000';
Ala.Parametro.PF      = 'trap';
Ala.Parametro.tor_d   = 0;
Ala.Parametro.theta0  = 0;
Ala.Parametro.theta0r = deg2rad(Ala.Parametro.theta0);
Ala.Parametro.tor     = @(y) Ala.Parametro.theta0r*(1-y^2/Ala.Parametro.b^2);

Ala.Parametro.Nc      = str2num(answer{10});
Ala.Parametro.nx      = Ala.Parametro.Nc+1;
Ala.Parametro.Nss     = str2num(answer{11});
Ala.Parametro.ny      = Ala.Parametro.Nss*2+1;
Ala.Parametro.Nt      = Ala.Parametro.Nc*Ala.Parametro.Nss*2;
Ala.Parametro.bias_x  = 1;
Ala.Parametro.bias_y  = 1;

CV.aoa     = transpose(cell2mat(textscan( answer{6}, '%f', 'Delimiter',',' )));
CV.aoar    = deg2rad(CV.aoa);
CV.beta   = str2num(answer{7});
CV.betar   = deg2rad(CV.beta);
CV.M       = str2num(answer{8});
CV.H       = str2num(answer{9});
[T,a,P,rho]= atmosisa(CV.H);
CV.U       = a*CV.M;
CV.q       = 0.5*rho*CV.U^2;
CV.naoa    = size(CV.aoar,2);


if strcmp(Ala.Parametro.PF,'trap') == 1
    C = @(y) Ctrap(y,Ala);

end

if str2double(Ala.Parametro.prof_R) < 10000
    Ala.Geometria.zc_r = @(x) NACA4(Ala.Parametro.prof_R,x);
elseif (str2double(Ala.Parametro.prof_R) >= 10^4) && (str2double(Ala.Parametro.prof_R) < 10^5)
    Ala.Geometria.zc_r = @(x) NACA5(Ala.Parametro.prof_R,x);
else
    error('NACA airfoil isn´t valid')
end

if str2double(Ala.Parametro.prof_T) < 10000
    Ala.Geometria.zc_t = @(x) NACA4(Ala.Parametro.prof_T,x);
elseif (str2double(Ala.Parametro.prof_T) >= 10^4) && (str2double(Ala.Parametro.prof_T) < 10^5)
    Ala.Geometria.zc_t = @(x) NACA5(Ala.Parametro.prof_T,x);
else
    error('NACA airfoil isn´t valid')
end

Ala.Geometria.C  = C;
Ala.Geometria.Cr = C(0);
Ala.Geometria.Ct = C(Ala.Parametro.b/2);

Ala.Geometria.Xca = 0.25*C(0)+2*tan(Ala.Parametro.Sweepr)/...
    Ala.Parametro.Sw*simpson(@(y)C(y)*y,0,Ala.Parametro.b/2,10);
Ala.Geometria.Yca = 2/Ala.Parametro.Sw*simpson(@(y)C(y)*y,0,Ala.Parametro.b/2,10);

Ala.Geometria.CMG = Ala.Parametro.Sw/Ala.Parametro.b;
Ala.Geometria.CMA = 2/Ala.Parametro.Sw*simpson(@(y)(C(y))^2,0,Ala.Parametro.b/2,10);


x0 = (Geom(0,1,Ala.Parametro.nx,Ala.Parametro.bias_x,true))';

y = (Geom(-Ala.Parametro.b/2,Ala.Parametro.b/2,...
    Ala.Parametro.ny,Ala.Parametro.bias_y,true))';
y(Ala.Parametro.Nss+1) = 0; 

Ala.Mesh.Y = ones(Ala.Parametro.nx,1)*y';
Ala.Mesh.X = zeros(Ala.Parametro.nx,Ala.Parametro.ny);
Ala.Mesh.Z = zeros(Ala.Parametro.nx,Ala.Parametro.ny);


for j = 1:Ala.Parametro.ny
    for i = 1:Ala.Parametro.nx
        yij = Ala.Mesh.Y(i,j);
        x0ij = x0(i);
        c = real(Ala.Geometria.C(yij));
        xij =Ala.Geometria.Cr/4+(x0ij-0.25)*c+abs(yij)*...
            tan(Ala.Parametro.Sweepr);
        zij = 0 ; 
        eps_y = Ala.Parametro.tor(yij);
        Ay = [cos(eps_y),sin(eps_y);-sin(eps_y),cos(eps_y)]; ...
            
        Ala.Mesh.X(i,j) = abs(yij)*tan(Ala.Parametro.Sweepr)+[1,0]*(Ay*...
            [xij-Ala.Geometria.Cr/4 - abs(yij)*...
            tan(Ala.Parametro.Sweepr),zij]')+Ala.Geometria.Cr/4;
        Ala.Mesh.Z(i,j) = [0,1]*(Ay*[xij-Ala.Geometria.Cr/4 - abs(yij)*...
            tan(Ala.Parametro.Sweepr),zij]')+abs(yij)*...
            tan(Ala.Parametro.dihedral); 
    end
end



Xs = zeros(Ala.Parametro.nx,Ala.Parametro.ny);
Zs = zeros(Ala.Parametro.nx,Ala.Parametro.ny);
for j = 1:Ala.Parametro.ny
    for i = 1:Ala.Parametro.nx
        yij = Ala.Mesh.Y(i,j);
        x0ij = x0(i);
        c = real(Ala.Geometria.C(yij));
        xij =Ala.Geometria.Cr/4+(x0ij-0.25)*c+abs(yij)*...
            tan(Ala.Parametro.Sweepr);
        zij = 0;
        eps_y = Ala.Parametro.tor(yij);
        Ay = [cos(eps_y),sin(eps_y);-sin(eps_y),cos(eps_y)]; ...
        Xs(i,j) = abs(yij)*tan(Ala.Parametro.Sweepr)+[1,0]*(Ay*[xij - ...
           Ala.Geometria.Cr/4 - abs(yij)*tan(Ala.Parametro.Sweepr),zij]')...
            +Ala.Geometria.Cr/4;
        Zs(i,j) = [0,1]*(Ay*[xij-Ala.Geometria.Cr/4 - abs(yij)*...
            tan(Ala.Parametro.Sweepr),zij]')+abs(yij)*...
            tan(Ala.Parametro.dihedralr);
        
       Ala.Geometria.torr(j) = eps_y;
       Ala.Geometria.tor(j)  = eps_y*180/pi;
    end
end

Xc     = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
Yc     = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
Zc     = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
dZc    = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
deltaP = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
Nvec   = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss,3);

for j = 1:2*Ala.Parametro.Nss
    for i = 1:Ala.Parametro.Nc
        
        Yc(i,j) = 0.5*(Ala.Mesh.Y(i,j)+Ala.Mesh.Y(i,j+1));
        Xc(i,j) = 0.5*(Xs(i,j)+Xs(i,j+1))+3/8*(Xs(i+1,j)-Xs(i,j)+...
            Xs(i+1,j+1)-Xs(i,j+1));
        Zc(i,j) = 0.5*(Zs(i,j)+Zs(i,j+1))+3/8*(Zs(i+1,j)-...
            Zs(i,j)+Zs(i+1,j+1)-Zs(i,j+1));
        
        epsilon = Ala.Geometria.C(Yc(i,j))*0.05; ...
            
        DeltaX = 2*epsilon;
        
        x1 = (Xc(i,j)-epsilon-tan(Ala.Parametro.Sweepr)*abs(Yc(i,j))...
            - Ala.Geometria.C(0)/4)/Ala.Geometria.C(Yc(i,j))+0.25;
        
        x2 = (Xc(i,j)+epsilon-tan(Ala.Parametro.Sweepr)*abs(Yc(i,j))-...
            Ala.Geometria.C(0)/4)/Ala.Geometria.C(Yc(i,j))+0.25;
        
        z1 = interp1([0,Ala.Parametro.b/2],...
            [Ala.Geometria.zc_r(x1),Ala.Geometria.zc_t(x1)],...
            abs(Yc(i,j)))*Ala.Geometria.C(Yc(i,j));
        
        z2 = interp1([0,Ala.Parametro.b/2],...
            [Ala.Geometria.zc_r(x2),Ala.Geometria.zc_t(x2)],...
            abs(Yc(i,j)))*Ala.Geometria.C(Yc(i,j));
        
        DeltaZ = z2-z1;
        dZc(i,j) = DeltaZ/DeltaX;
       
        deltaP(i,j) = atan(dZc(i,j));
        
        v1 = [Xc(i,j),Yc(i,j),Zc(i,j)]-[Xs(i,j),Ala.Mesh.Y(i,j),Zs(i,j)];
        v2 = [Xc(i,j),Yc(i,j),Zc(i,j)]-[Xs(i,j+1),Ala.Mesh.Y(i,j+1),Zs(i,j+1)];
        nvec = cross(v2,v1)/norm(cross(v2,v1));
        
        A = [cos(deltaP(i,j)),0,-sin(deltaP(i,j));...
            0,1,0;...
            sin(deltaP(i,j)),0,cos(deltaP(i,j))];
        Nvec(i,j,:) = (A*nvec')';
    end
end





Xt = zeros(Ala.Parametro.Nc,Ala.Parametro.ny);
Zt = zeros(Ala.Parametro.Nc,Ala.Parametro.ny);

for i = 1:Ala.Parametro.Nc
    for j = 1:Ala.Parametro.ny
        Xt(i,j) = 0.75*Xs(i,j)+0.25*Xs(i+1,j);
        Zt(i,j) = 0.75*Zs(i,j)+0.25*Zs(i+1,j);
    end
end



figure
hold on
mesh(Xs,Ala.Mesh.Y,Zs)
plot3(Xc,Yc,Zc,'.k')
quiver3(Xc,Yc,Zc,Nvec(:,:,1),Nvec(:,:,2),Nvec(:,:,3),0.5)
plot3(Xt',(Ala.Mesh.Y(1:Ala.Parametro.Nc,:))',Zt','--r')
hold off
axis('equal')
grid


Xt = [Xt; Xt(Ala.Parametro.Nc,:) + (Xt(Ala.Parametro.Nc,:) -  ...
    Xt(Ala.Parametro.Nc-1,:))*1000] ;
Zt = [Zt; Zt(Ala.Parametro.Nc,:) + (Zt(Ala.Parametro.Nc,:) - ...
    Zt(Ala.Parametro.Nc-1,:))*1000] ;


Ala.Mesh.X         = Ala.Mesh.X + Ala.Parametro.xoffset;
Ala.Mesh.Y         = Ala.Mesh.Y + Ala.Parametro.yoffset;
Ala.Mesh.Z         = Ala.Mesh.Z + Ala.Parametro.zoffset;
Ala.Mesh.Control.X = Xc + Ala.Parametro.xoffset;
Ala.Mesh.Control.Y = Yc + Ala.Parametro.yoffset;
Ala.Mesh.Control.Z = Zc + Ala.Parametro.zoffset;
Ala.Mesh.Node.X    = Xt + Ala.Parametro.xoffset;
Ala.Mesh.Node.Y    = Ala.Mesh.Y;
Ala.Mesh.Node.Z    = Zt + Ala.Parametro.zoffset;

Ala.Mesh.Xs        = Xs + Ala.Parametro.xoffset;
Ala.Mesh.Zs        = Zs + Ala.Parametro.zoffset;


Ala.Mesh.Nvec      = Nvec;
XC    = Ala.Mesh.Control.X;
YC    = Ala.Mesh.Control.Y;
ZC    = Ala.Mesh.Control.Z;


Xn    = Ala.Mesh.Node.X;
Yn    = Ala.Mesh.Node.Y;
Zn    = Ala.Mesh.Node.Z;



Xs    = Ala.Mesh.Xs;
Zs    = Ala.Mesh.Zs;
Nvec  = Ala.Mesh.Nvec;



Aij = zeros(Ala.Parametro.Nt,Ala.Parametro.Nt);


for i1 = 1:Ala.Parametro.Nc
    for j1 = 1:2*Ala.Parametro.Nss
                    PC = [XC(i1,j1), YC(i1,j1), ZC(i1,j1)];
        for i2 = 1:Ala.Parametro.Nc
            for j2 = 1:2*Ala.Parametro.Nss
                
                    
                    P1 = [Xn(i2,j2), Yn(i2,j2), Zn(i2,j2)];
                    P2 = [Xn(i2,j2+1), Yn(i2,j2+1), Zn(i2,j2+1)];
                    P3 = [Xn(i2+1,j2+1), Yn(i2+1,j2+1), Zn(i2+1,j2+1)];
                    P4 = [Xn(i2+1,j2), Yn(i2+1,j2), Zn(i2+1,j2)];

                    
                    
                
                
                r21 = P2 - P1;
                r1  = PC - P1;
                r2  = PC - P2;
                
           
                psi   = cross(r1,r2)/((norm(cross(r1,r2)))^2); 
                omega = r21*(r1/norm(r1)-r2/norm(r2))'; 
                V12   = 1/(4*pi)*psi*omega;
                
                if norm(cross(r1/norm(r1),r2/norm(r2))) < 1e-12
                    V12 = 0;
                end
                
                
                
                r32 = P3 - P2;
                r3 = PC - P3;

                
                psi = cross(r2,r3)/((norm(cross(r2,r3)))^2); 
                omega = r32*(r2/norm(r2)-r3/norm(r3))';
                V23 = 1/(4*pi)*psi*omega;
                
                 if norm(cross(r2/norm(r2),r3/norm(r3))) < 1e-12
                    V23 = 0;
                 end
                
                 
                
                
                
                r43 = P4 - P3;
                r4 = PC - P4;

                
                psi = cross(r3,r4)/((norm(cross(r3,r4)))^2); 
                omega = r43*(r3/norm(r3)-r4/norm(r4))'; 
                V34 = 1/(4*pi)*psi*omega;
                
                
                if norm(cross(r3/norm(r3),r4/norm(r4))) < 1e-12
                    V34 = 0;
                end
                
                
                
                
                
                r14 = P1 - P4;

                psi = cross(r4,r1)/((norm(cross(r4,r1)))^2); 
                omega = r14*(r4/norm(r4)-r1/norm(r1))'; 
                V41 = 1/(4*pi)*psi*omega;
                if norm(cross(r4/norm(r4),r1/norm(r1))) < 1e-12
                    V41 = 0;
                end
                
                
                
                Vind = V12+V23+V34+V41;
                
              
                
                vn = Vind*[Nvec(i1,j1,1);Nvec(i1,j1,2);Nvec(i1,j1,3)];

                i = i1+Ala.Parametro.Nc*(j1-1);
                j = i2+Ala.Parametro.Nc*(j2-1);
                Aij(i,j) = Vind(end);


            end
        end
    end
end

Ala.VML.AIC     = Aij;




for k = 1:CV.naoa
    
    
    if CV.aoa(k) <0
        
        CV.aoa(k)         = -CV.aoa(k);
        name              = sprintf('AoAneg%1i', 100*CV.aoa(k));
        name              = strrep(name,'.','');
        CV.aoa(k)         = -CV.aoa(k);
    else
        
        name              = sprintf('AoA%1i', 100*CV.aoa(k));
        name              = strrep(name,'.','');
    end
        
    Ala.VML.(name).clj   = zeros(1,2*Ala.Parametro.Nss);
    
    Talpha = [cos(CV.aoar(k)),0,-sin(CV.aoar(k));...
        0,1,0;...
        sin(CV.aoar(k)),0,cos(CV.aoar(k))];
    Tbeta  = [cos(CV.betar),sin(CV.betar),0;...
        -sin(CV.betar),cos(CV.betar),0;...
        0,0,1];
    
    U      = (CV.U*Tbeta*Talpha*[1,0,0]')';
    
    Ala.VML.(name).Wyij   = zeros(Ala.Parametro.Nc*2*Ala.Parametro.Nss,1);
    
    for i1 = 1:Ala.Parametro.Nc
        for j1 = 1:2*Ala.Parametro.Nss
            vinf_n = U*[Nvec(i1,j1,1);Nvec(i1,j1,2);Nvec(i1,j1,3)];
            i = i1+Ala.Parametro.Nc*(j1-1);
            Ala.VML.(name).Wyij(i) = -vinf_n;
        end
    end
    

    
    Ala.VML.(name).Gamma = Aij\Ala.VML.(name).Wyij;
    
    
    
    figures = true;
    

    
    Gamma2 = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
    for j = 1:2*Ala.Parametro.Nss
        for i = 1:Ala.Parametro.Nc
            Gamma2(i,j) = Ala.VML.(name).Gamma(i+Ala.Parametro.Nc*(j-1));
        end
    end
    
    Gamma2 = [Gamma2(1,:) ; Gamma2(2:end,:) - Gamma2(1:end-1,:)] ;
    
    if figures == true
        esc = 1/max(max(abs(Gamma2))); 
        figure
        
        
        hold on
        surf(Xs,Ala.Mesh.Y,Zs,Gamma2)
        
        colorbar
        if Ala.Parametro.Nc == 1
            mesh ([Ala.Mesh.Control.X;Ala.Mesh.Control.X],...
                [Ala.Mesh.Control.Y;Ala.Mesh.Control.Y],...
                esc*[Gamma2;Gamma2],[Gamma2;Gamma2])
        else
            mesh (Ala.Mesh.Control.X,Ala.Mesh.Control.Y,esc*Gamma2,Gamma2,'facecol','none')
        end
        hold off
        az = 90;
        el = 90;
        view(az, el);
       
        xlabel('Chordwise [m]')
        ylabel('Spanwise [m]')
        title('$\Gamma$')
        grid
        set(gca, 'DataAspectRatio', [1, 1, 1/Ala.Parametro.b])
    end
    
       
    Ala.VML.(name).clij = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
    Ala.VML.(name).cdij = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
    for i = 1:Ala.Parametro.Nc
        for j = 1:2*Ala.Parametro.Nss
            n = [Nvec(i,j,1),Nvec(i,j,2),Nvec(i,j,3)];
            cij = 0.5*(Xs(i+1,j)-Xs(i,j)+Xs(i+1,j+1)-Xs(i,j+1));
            Ala.VML.(name).Cp = 2*Gamma2(i,j)/CV.U;
            Ala.VML.(name).cdij(i,j) = Ala.VML.(name).Cp/cij*U/CV.U*n';
            if Ala.VML.(name).Cp == 0
                Ala.VML.(name).clij(i,j) = 0;
            else
                Ala.VML.(name).clij(i,j) = Ala.VML.(name).Cp/...
                    abs(Ala.VML.(name).Cp)*norm(Ala.VML.(name).Cp/cij*n-...
                    Ala.VML.(name).cdij(i,j)*U/CV.U);
            end
        end
    end
    
Ala.VML.(name).Cpij = 2*Gamma2/CV.U;
   
    if figures == true
        esc2 = 1/max(max(abs(Ala.VML.(name).clij)));
        figure
        
        hold on
        surf(Ala.Mesh.X,Ala.Mesh.Y,Ala.Mesh.Z,Ala.VML.(name).clij)
    
        colorbar
        if Ala.Parametro.Nc == 1
            mesh ([Ala.Mesh.Control.X;Ala.Mesh.Control.X],...
                [Ala.Mesh.Control.Y;Ala.Mesh.Control.Y],...
                esc2*[Ala.VML.(name).clij;Ala.VML.(name).clij],...
                [Ala.VML.(name).clij;Ala.VML.(name).clij])
        else
            mesh(Ala.Mesh.Control.X,Ala.Mesh.Control.Y,...
                esc2*Ala.VML.(name).clij,Ala.VML.(name).clij,'facecol','none')
        end
        hold off
        az = 90;
        el = 90;
        view(az, el);
      
        xlabel('Chordwise [m]')
        ylabel('Spanwise [m]')
        title('$C_{lij}$')
        grid
        set(gca, 'DataAspectRatio', [1, 1, 1/Ala.Parametro.b])
    end
    
    Z0     = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
    
    if figures == true
        esc3 = 1/max(max(abs(Ala.VML.(name).Cpij)));
        figure

        hold on
        surf(Ala.Mesh.X,Ala.Mesh.Y,Ala.Mesh.Z,Ala.VML.(name).Cpij)
        colorbar
        if Ala.Parametro.Nc == 1
            mesh ([Ala.Mesh.Control.X;Ala.Mesh.Control.X],...
                [Ala.Mesh.Control.Y;Ala.Mesh.Control.Y],...
                esc3*[Ala.VML.(name).Cpij;Ala.VML.(name).Cpij],...
                [Ala.VML.(name).Cpij;Ala.VML.(name).Cpij])
        else
            mesh(Ala.Mesh.Control.X,Ala.Mesh.Control.Y,...
                esc3*Ala.VML.(name).Cpij,Ala.VML.(name).Cpij,'facecol','none')
        end
        hold off
        az = 90;
        el = 90;
        view(az, el);
    
        xlabel('Chordwise [m]')
        ylabel('Spanwise [m]')
        title('$C_{p}$')
        grid
        set(gca, 'DataAspectRatio', [1, 1, 1/Ala.Parametro.b])
    end
    
figure()
surf(Ala.Mesh.Control.X,Ala.Mesh.Control.Y,Gamma2); hold on
mesh(Ala.Mesh.Control.X,Ala.Mesh.Control.Y,Z0,'facecol','none')
colorbar

az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$\Gamma$')
grid on

figure()
surf(Ala.Mesh.Control.X,Ala.Mesh.Control.Y,Ala.VML.(name).Cpij); hold on
mesh(Ala.Mesh.Control.X,Ala.Mesh.Control.Y,Z0,'facecol','none')
colorbar

az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$C_{p}$')
grid on

figure()
surf(Ala.Mesh.Control.X(:,Ala.Parametro.Nss+1:end),Ala.Mesh.Control.Y(:,Ala.Parametro.Nss+1:end),Ala.VML.(name).Cpij(:,Ala.Parametro.Nss+1:end)); hold on
mesh(Ala.Mesh.Control.X(:,Ala.Parametro.Nss+1:end),Ala.Mesh.Control.Y(:,Ala.Parametro.Nss+1:end),Z0(:,Ala.Parametro.Nss+1:end),'facecol','none')
colorbar

az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$C_{p}$')
grid on    
    
    
   
    
    clij2 = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
    cd    = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
    
    for j = 1:2*Ala.Parametro.Nss
        for i = 1:Ala.Parametro.Nc
            cij = 0.5*(Xs(i+1,j)-Xs(i,j)+Xs(i+1,j+1)-Xs(i,j+1));
            clij2(i,j) = Ala.VML.(name).clij(i,j)*cij/...
                Ala.Geometria.C(Ala.Mesh.Control.Y(i,j));
            cd(i,j) = Ala.VML.(name).cdij(i,j)*cij/...
                Ala.Geometria.C(Ala.Mesh.Control.Y(i,j));
        end
        Ala.VML.(name).clj(j) = sum(clij2(:,j));
        Ala.VML.(name).cdy(j) = sum(cd(:,j));
    end



if figures == true
        figure

        plot([-Ala.Parametro.b/2,Ala.Mesh.Control.Y(1,:),Ala.Parametro.b/2],...
            [0,Ala.VML.(name).clj,0])
        xlabel('Spanwise (m)')
        ylabel('$c_l(y)$')
        title('$C_{l}$')
        grid
end



figure()
plot(Ala.Mesh.Control.Y(1,Ala.Parametro.Nss+1:end),Ala.VML.(name).clj(Ala.Parametro.Nss+1:end))
xlabel('Spanwise (m)')
ylabel('$c_l(y)$')
title('$C_{l}$')
grid on



c14   = floor(Ala.Parametro.Nc/4);
c34   = round(3*Ala.Parametro.Nc/4);
figure()
plot(Ala.Mesh.Control.Y(1,Ala.Parametro.Nss+1:end),Ala.VML.(name).Cpij(1,Ala.Parametro.Nss+1:end)); hold on
plot(Ala.Mesh.Control.Y(c14,Ala.Parametro.Nss+1:end),Ala.VML.(name).Cpij(c14,Ala.Parametro.Nss+1:end))
plot(Ala.Mesh.Control.Y(c34,Ala.Parametro.Nss+1:end),Ala.VML.(name).Cpij(c34,Ala.Parametro.Nss+1:end));
xlabel('Spanwise (m)')
ylabel('$c_p$')
leg=legend('Leading Edge','Chord 1/4','Chord 3/4', 'Trailing Edge');
set(leg,'interpreter','latex','location','best')
grid on
title('$C_{p}$')

s14   = floor(Ala.Parametro.Nss+Ala.Parametro.Nss/4+1);
s34   = round(Ala.Parametro.Nss+3*Ala.Parametro.Nc/4+1);
figure()
plot(Ala.Mesh.Control.X(:,1),Ala.VML.(name).Cpij(:,1)); hold on
plot(Ala.Mesh.Control.X(:,s14),Ala.VML.(name).Cpij(:,s14))
plot(Ala.Mesh.Control.X(:,s34),Ala.VML.(name).Cpij(:,s34));
xlabel('Chordwise (m)')
ylabel('$c_p$')
leg=legend('Root','1/4 Semi Span','3/4 Semi Span', 'Trailing Edge');
set(leg,'interpreter','latex','location','best')
grid on
title('$C_{p}$')





    
bj = zeros(1,2*Ala.Parametro.Nss);
bcj = zeros(1,2*Ala.Parametro.Nss);
for j = 1:2*Ala.Parametro.Nss
        bj(j) = Ala.Mesh.Y(1,j+1)-Ala.Mesh.Y(i,j);
        bcj(j) = bj(j)*Ala.Geometria.C(YC(1,j));
end
Ala.VML.(name).cL  = 1/Ala.Parametro.Sw*Ala.VML.(name).clj*bcj';
Ala.VML.(name).cDi = 1/Ala.Parametro.Sw*Ala.VML.(name).cdy*bcj';
Ala.VML.(name).E   = Ala.VML.(name).cL/Ala.VML.(name).cDi;



CLplot(k)           = Ala.VML.(name).cL;
CDplot(k)           = Ala.VML.(name).cDi;
Eplot(k)            = Ala.VML.(name).E;


    
m0ij = zeros(Ala.Parametro.Nc,2*Ala.Parametro.Nss);
    
for j = 1:2*Ala.Parametro.Nss
   for i = 1:Ala.Parametro.Nc
            cij = 0.5*(Xs(i+1,j)-Xs(i,j)+Xs(i+1,j+1)-Xs(i,j+1));
            Xc4ij = 0.5*(Xn(i,j)+Xn(i,j+1));
            m0ij(i,j) = -Ala.VML.(name).clij(i,j)*Xc4ij*cij*bj(j);
   end
end
Ala.VML.(name).cM0y = 1/(Ala.Parametro.Sw*Ala.Geometria.C(0))*sum(sum(m0ij));


end


if CV.naoa > 1
    
    Cla  = CLplot./CV.aoar;
    
    figure()
    plot(CV.aoa,CLplot)
    xlabel('$\alpha$')
    ylabel('$C_L$')
    grid on
    title('$C_{L}$ vs $\alpha$')
    
    figure()
    plot(CV.aoa,CDplot)
    xlabel('$\alpha$')
    ylabel('$C_D$')
    grid on
    title('$C_{D}$ vs $\alpha$')
   
    figure()
    plot(CLplot,CDplot)
    xlabel('$C_L$')
    ylabel('$C_D$')
    grid on
    title('Polar')
    
end
