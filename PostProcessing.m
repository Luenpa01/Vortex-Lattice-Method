function PostProcessing(Ala,Tail,CV)

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
        surf(Ala.Mesh.Xs,Ala.Mesh.Y,Ala.Mesh.Zs,Gamma2)
        surf(Tail.Mesh.Xs,Tail.Mesh.Y,Tail.Mesh.Zs,Gamma2)

        colorbar
        if Ala.Parametro.Nc == 1
            mesh ([Ala.Mesh.Control.X;Ala.Mesh.Control.X],...
                [Ala.Mesh.Control.Y;Ala.Mesh.Control.Y],...
                esc*[Gamma2;Gamma2],[Gamma2;Gamma2])
            mesh ([Tail.Mesh.Control.X;Tail.Mesh.Control.X],...
                [Tail.Mesh.Control.Y;Tail.Mesh.Control.Y],...
                esc*[Gamma2;Gamma2],[Gamma2;Gamma2])
        else
            mesh (Ala.Mesh.Control.X,Ala.Mesh.Control.Y,esc*Gamma2,Gamma2,'facecol','none')
            mesh (Tail.Mesh.Control.X,Tail.Mesh.Control.Y,esc*Gamma2,Gamma2,'facecol','none')
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
    
     
    if figures == true
        esc2 = 1/max(max(abs(Ala.VML.(name).clij)));
        figure
        
        hold on
        surf(Ala.Mesh.X,Ala.Mesh.Y,Ala.Mesh.Z,Ala.VML.(name).clij)
        surf(Tail.Mesh.X,Tail.Mesh.Y,Tail.Mesh.Z,Tail.VLM.(name).clij)
   
        colorbar
        if Ala.Parametro.Nc == 1
            mesh ([Ala.Mesh.Control.X;Ala.Mesh.Control.X],...
                [Ala.Mesh.Control.Y;Ala.Mesh.Control.Y],...
                esc2*[Ala.VML.(name).clij;Ala.VML.(name).clij],...
                [Ala.VML.(name).clij;Ala.VML.(name).clij])
            mesh ([Tail.Mesh.Control.X;Tail.Mesh.Control.X],...
                [Tail.Mesh.Control.Y;Tail.Mesh.Control.Y],...
                esc2*[Tail.VLM.(name).clij;Tail.VLM.(name).clij],...
                [Tail.VLM.(name).clij;Tail.VLM.(name).clij])
        else
            mesh(Ala.Mesh.Control.X,Ala.Mesh.Control.Y,...
                esc2*Ala.VML.(name).clij,Ala.VML.(name).clij,'facecol','none')
            mesh(Tail.Mesh.Control.X,Tail.Mesh.Control.Y,...
                esc2*Tail.VLM.(name).clij,Tail.VLM.(name).clij,'facecol','none')
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
        surf(Tail.Mesh.X,Tail.Mesh.Y,Tail.Mesh.Z,Tail.VLM.(name).Cpij)

        colorbar
        if Ala.Parametro.Nc == 1
            mesh ([Ala.Mesh.Control.X;Ala.Mesh.Control.X],...
                [Ala.Mesh.Control.Y;Ala.Mesh.Control.Y],...
                esc3*[Ala.VML.(name).Cpij;Ala.VML.(name).Cpij],...
                [Ala.VML.(name).Cpij;Ala.VML.(name).Cpij])
            mesh ([Tail.Mesh.Control.X;Tail.Mesh.Control.X],...
                [Tail.Mesh.Control.Y;Tail.Mesh.Control.Y],...
                esc3*[Tail.VLM.(name).Cpij;Tail.VLM.(name).Cpij],...
                [Tail.VLM.(name).Cpij;Tail.VLM.(name).Cpij])
        else
            mesh(Ala.Mesh.Control.X,Ala.Mesh.Control.Y,...
                esc3*Ala.VML.(name).Cpij,Ala.VML.(name).Cpij,'facecol','none')
            mesh(Tail.Mesh.Control.X,Tail.Mesh.Control.Y,...
                esc3*Tail.VLM.(name).Cpij,Tail.VLM.(name).Cpij,'facecol','none')
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
hold off

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
hold off

figure()
surf(Ala.Mesh.Control.X(:,Ala.Parametro.Nss+1:end),Ala.Mesh.Control.Y(:,Ala.Parametro.Nss+1:end),Ala.VML.(name).Cpij(:,Ala.Parametro.Nss+1:end)); hold on
surf(Tail.Mesh.Control.X(:,Tail.Parameters.Nss+1:end),Tail.Mesh.Control.Y(:,Tail.Parameters.Nss+1:end),Tail.VLM.(name).Cpij(:,Tail.Parameters.Nss+1:end))
mesh(Ala.Mesh.Control.X(:,Ala.Parametro.Nss+1:end),Ala.Mesh.Control.Y(:,Ala.Parametro.Nss+1:end),Z0(:,Ala.Parametro.Nss+1:end),'facecol','none')
mesh(Tail.Mesh.Control.X(:,Tail.Parameters.Nss+1:end),Tail.Mesh.Control.Y(:,Tail.Parameters.Nss+1:end),Z0(:,Tail.Parameters.Nss+1:end),'facecol','none')
colorbar

az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$C_{p}$')
grid on    
    

end

end