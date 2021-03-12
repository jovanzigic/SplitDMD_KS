function [ROM_error_DMD,ROM_error_OD] = sim_DMD_KS(realXdmd,errorDMD,realXod,errorOD,r,interval,M)

    Nx = size(realXdmd,1);
    Nt = size(realXdmd,2);

    x = linspace(0,1,Nx);
    t = linspace(0,interval,Nt);
    
    % Simulate DMD reduced model
    figure
    mesh(x,t,realXdmd'); 
    %view([.7 -.8 .6])
    view([0 0 1])
    title(['Standard DMD Model, order = ' , num2str(r) ] )
    xlabel('x'); ylabel('t'); zlabel('z_{DMD}')
    xlim([0 1]) 
    ylim([0 interval])
%     zlim([-20 20])
    
    % Simulate Error of DMD ROM
    figure
    mesh(x,t,errorDMD')
    %view([.7 -.8 .6])
    view([0 0 1])
    xlabel('x'); ylabel('t'); zlabel('z-z_{DMD}')
    title(['Error in Standard DMD Model, order = ' , num2str(r)] )
    xlim([0 1]) 
    ylim([0 interval])
%     zlim([-20 20])
    
    % Simulate Optimized DMD reduced model
    figure
    mesh(x,t,realXod'); 
    %view([.7 -.8 .6])
    view([0 0 1])
    title(['Optimized DMD Model, order = ' , num2str(r) ] )
    xlabel('x'); ylabel('t'); zlabel('z_{OD}')
    xlim([0 1]) 
    ylim([0 interval])
%     zlim([-20 20])
    
    % Simulate Error of Optimized DMD ROM
    figure
    mesh(x,t,errorOD')
    %view([.7 -.8 .6])
    view([0 0 1])
    xlabel('x'); ylabel('t'); zlabel('z-z_{OD}')
    title(['Error in Optimized DMD Model, order = ' , num2str(r)] )
    xlim([0 1]) 
    ylim([0 interval])
%     zlim([-20 20])

    % L2 Norm 

    delta_t = (t(2)-t(1))/Nt;

    ROM_error_DMD = 0.5*errorDMD(:,1)'*M*errorDMD(:,1);
    for i=2:Nt-1
    ROM_error_DMD = ROM_error_DMD + errorDMD(:,i)'*M*errorDMD(:,i);
    end
    ROM_error_DMD = ROM_error_DMD + 0.5*errorDMD(:,Nt)'*M*errorDMD(:,Nt);

    ROM_error_DMD = sqrt(delta_t*ROM_error_DMD); 

    ROM_error_OD = 0.5*errorOD(:,1)'*M*errorOD(:,1);
    for i=2:Nt-1
    ROM_error_OD = ROM_error_OD + errorOD(:,i)'*M*errorOD(:,i);
    end
    ROM_error_OD = ROM_error_OD + 0.5*errorOD(:,Nt)'*M*errorOD(:,Nt);

    ROM_error_OD = sqrt(delta_t*ROM_error_OD); 

end