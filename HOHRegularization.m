clc; clearvars; 

% Dados fornecidos Ghia
% coordenadas
y_ghia = [1.0, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, ...
     0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0];
x_ghia = [1.0, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, ...
      0.5, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0];
% velocidades para Re = 10000
u_ghia = [1.0, 0.47221, 0.47783, 0.48070, 0.47804, 0.34635, 0.20673, 0.08344, ...
    0.03111, -0.07540, -0.23186, -0.32709, -0.38000, -0.41657, ...
     -0.42537, -0.42735, 0.0];
v_ghia = [0.0, -0.54302, -0.52987, -0.49099, -0.45863, -0.41496, -0.36737, ...
     -0.30719, 0.00831, 0.27224, 0.28003, 0.35070, 0.41487, ...
     0.43124, 0.43733, 0.43983, 0.0];
% velocidades para Re = 100
% v_ghia = [0.0, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, ...
%    -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, ...
%    0.10091, 0.09233, 0.0]; 
% u_ghia = [1.0, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, ...
%         0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, ...
%        -0.06434, -0.04775, -0.04192, -0.03717, 0.0];

% GENERAL FLOW CONSTANTS 
lx = 129; 
ly = lx; 

uLid  = .1/sqrt(3); % horizontal lid velocity 
vLid  = 0;    % vertical lid velocity 
Re    = 10000;  % Reynolds number 100 or 10000
nu    = uLid *(lx-1) / Re;     % kinematic viscosity 
omega = 1. / (3*nu+1./2.); % relaxation parameter 
maxT  = 1000000000; % total number of iterations 500000
tPlot = 4000;    % cycles for graphical output4000

% D2Q9 LATTICE CONSTANTS 
w   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36]; 
cx  = [ 0, 1, 0, -1, 0, 1, -1, -1, 1]; 
cy  = [ 0, 0, 1, 0, -1, 1, 1, -1, -1]; 
opp = [ 1, 4, 5, 2,  3, 8, 9,  6,  7]; 
lid = [2: (lx-1)]; 

% INITIAL CONDITION: (rho=1, u=0) ==> fIn(i) = t(i) 
fIn = reshape( w' * ones(1,lx*ly), 9, lx, ly); 

% MAIN LOOP (TIME CYCLES) 
tic

%video
%videoNome = 'HORR.mp4'; 
%vidObj = VideoWriter(videoNome, 'MPEG-4');
%vidObj.FrameRate = 30; % utlizar 30 ou 60
%open(vidObj);
%hFig = figure;

for cycle = 1:maxT 

% MACROSCOPIC VARIABLES 
rho = sum(fIn);
ux = reshape ( (cx * reshape(fIn,9,lx*ly)), 1,lx,ly ) ./rho; 
uy = reshape ( (cy * reshape(fIn,9,lx*ly)), 1,lx,ly ) ./rho; 

% MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS 
ux(:,lid,ly) = uLid; %lid x - velocity 
uy(:,lid,ly) = vLid; %lid y - velocity 

% BOUNDARY VELOCITY CONDITION - BOTTOM WALL
ux(:,:,1) = 0.; 
uy(:,:,1) = 0.; 

% BOUNDARY VELOCITY CONDITION - TOP LEFT WALL
ux(:,1,ly) = 0.; 
uy(:,1,ly) = 0.; 

% BOUNDARY VELOCITY CONDITION - TOP RIGHT WALL
ux(:,lx,ly) = 0.; 
uy(:,lx,ly) = 0.; 

% BOUNDARY VELOCITY CONDITION - LEFT WALL
ux(:,1,:) = 0.; 
uy(:,1,:) = 0.; 

% BOUNDARY VELOCITY CONDITION - RIGHT WALL
ux(:,lx,:) = 0.; 
uy(:,lx,:) = 0.; 

% REGULARIZED BOUNDARY CONDITIONS - BOTTOM WALL 
rhos = zeros(1,lx-2);
rhos_m20Istar = zeros(1,lx-2);
rhos_m11Istar = zeros(1,lx-2);
rhos_m02Istar = zeros(1,lx-2);

for i=1:9
   if cy(i)<=0 
      rhos = rhos + fIn(i,2:(lx-1),1);      
      rhos_m20Istar = rhos_m20Istar + fIn(i,2:(lx-1),1)*(cx(i)*cx(i)-1/3); 
      rhos_m11Istar = rhos_m11Istar + fIn(i,2:(lx-1),1)*(cx(i)*cy(i)); 
      rhos_m02Istar = rhos_m02Istar + fIn(i,2:(lx-1),1)*(cy(i)*cy(i)-1/3);
   end
end

rhoprime = (12*rhos+(1-omega)*9*rhos_m02Istar)/(9+omega);
rho_m20prime = 6/5*rhos_m20Istar;
rho_m11prime = 2*rhos_m11Istar;
rho_m02prime = (2*rhos+15*rhos_m02Istar)/(9+omega);

for i=1:9 
    fIn(i,2:(lx-1),1) =  w(i) * ... 
        ( rhoprime + ...
        9/2*(cx(i)*cx(i)-1/3)*rho_m20prime + ...
        9/2*(cy(i)*cy(i)-1/3)*rho_m02prime + ...    
        9*cx(i)*cy(i)*rho_m11prime); 
       
end

% END OF BOUNDARY CONDITIONS - BOTTOM WALL
    
% REGULARIZED BOUNDARY CONDITIONS - TOP WALL 
rhos = zeros(1,lx-2);
rhos_m20Istar = zeros(1,lx-2);
rhos_m11Istar = zeros(1,lx-2);
rhos_m02Istar = zeros(1,lx-2);

for i=1:9
   if cy(i)>=0 
      rhos = rhos + fIn(i,2:(lx-1),ly);      
      rhos_m20Istar = rhos_m20Istar + fIn(i,2:(lx-1),ly)*(cx(i)*cx(i)-1/3); 
      rhos_m11Istar = rhos_m11Istar + fIn(i,2:(lx-1),ly)*(cx(i)*cy(i)); 
      rhos_m02Istar = rhos_m02Istar + fIn(i,2:(lx-1),ly)*(cy(i)*cy(i)-1/3);
   end
end

rhoprime = (12*rhos+(1-omega)*9*rhos_m02Istar)/(9+omega);
rho_m20prime = (12*rhos_m20Istar - ...
    24*rhos_m11Istar*uLid + ...
    7*rhoprime*uLid*uLid + ...
    27 * rhos_m02Istar*(uLid^2))/10.0;
rho_m11prime = (1/2)*(4*rhos_m11Istar-rhoprime*uLid - 3.0*rhos_m02Istar*(uLid));
rho_m02prime = (1/6)*(rhoprime +(9*rhos_m02Istar));

rhoprime_ux = rhoprime * uLid;

for i=1:9 
    fIn(i,2:(lx-1),ly) =  w(i) * ... 
        ( rhoprime + 3*cx(i)*rhoprime_ux +...
        9/2*(cx(i)*cx(i)-1/3)*rho_m20prime + ...
        9/2*(cy(i)*cy(i)-1/3)*rho_m02prime + ...    
        9*cx(i)*cy(i)*rho_m11prime + ...
        27/2 * (2 * uLid * rho_m11prime) * (cx(i)^2 - 1/3) * cy(i) + ...
        27/2 * (uLid * rho_m02prime) * (cy(i)^2 - 1/3) * cx(i));
end

% END OF BOUNDARY CONDITIONS - TOP WALL
    
% REGULARIZED BOUNDARY CONDITIONS - BOTTOM LEFT WALL 
rhos = 0;
rhos_m20Istar = 0;
rhos_m11Istar = 0;
rhos_m02Istar = 0;

for i=1:9
   if (cx(i)<=0 && cy(i)<=0)
      rhos = rhos + fIn(i,1,1);      
      rhos_m20Istar = rhos_m20Istar + fIn(i,1,1)*(cx(i)*cx(i)-1/3); 
      rhos_m11Istar = rhos_m11Istar + fIn(i,1,1)*(cx(i)*cy(i)); 
      rhos_m02Istar = rhos_m02Istar + fIn(i,1,1)*(cy(i)*cy(i)-1/3);
   end
end

rhoprime = -12/(16+9*omega)*(-3*rhos-3*rhos_m02Istar+omega*3*rhos_m02Istar+...
            7*rhos_m11Istar-omega*7*rhos_m11Istar-3*rhos_m20Istar+...
            omega*3*rhos_m20Istar);
rho_m20prime = 2/9*(rhoprime-6*rhos_m11Istar+9*rhos_m20Istar);
rho_m11prime = 1/27*(-7*rhoprime-18*rhos_m02Istar+132*rhos_m11Istar-18*rhos_m20Istar);
rho_m02prime = 2/9*(rhoprime-6*rhos_m11Istar+9*rhos_m02Istar);

for i=1:9 
    fIn(i,1,1) =  w(i) * ... 
        ( rhoprime + ...
        9/2*(cx(i)*cx(i)-1/3)*rho_m20prime + ...
        9/2*(cy(i)*cy(i)-1/3)*rho_m02prime + ...    
        9*cx(i)*cy(i)*rho_m11prime);
end
% END OF BOUNDARY CONDITIONS - BOTTOM LEFT WALL

% REGULARIZED BOUNDARY CONDITIONS - BOTTOM RIGHT WALL 
rhos = 0;
rhos_m20Istar = 0;
rhos_m11Istar = 0;
rhos_m02Istar = 0;

for i=1:9
   if (cx(i)>=0 && cy(i)<=0)
      rhos = rhos + fIn(i,lx,1);      
      rhos_m20Istar = rhos_m20Istar + fIn(i,lx,1)*(cx(i)*cx(i)-1/3); 
      rhos_m11Istar = rhos_m11Istar + fIn(i,lx,1)*(cx(i)*cy(i)); 
      rhos_m02Istar = rhos_m02Istar + fIn(i,lx,1)*(cy(i)*cy(i)-1/3);
   end
end

rhoprime = -12/(16+9*omega)*(-3*rhos-3*rhos_m02Istar+omega*3*rhos_m02Istar-...
            7*rhos_m11Istar+omega*7*rhos_m11Istar-3*rhos_m20Istar+...
            omega*3*rhos_m20Istar);
rho_m20prime = 2/9*(rhoprime+6*rhos_m11Istar+9*rhos_m20Istar);
rho_m11prime = 1/27*(7*rhoprime+18*rhos_m02Istar+132*rhos_m11Istar+18*rhos_m20Istar);
rho_m02prime = 2/9*(rhoprime+6*rhos_m11Istar+9*rhos_m02Istar);

for i=1:9 
    fIn(i,lx,1) =  w(i) * ... 
        ( rhoprime + ...
        9/2*(cx(i)*cx(i)-1/3)*rho_m20prime + ...
        9/2*(cy(i)*cy(i)-1/3)*rho_m02prime + ...    
        9*cx(i)*cy(i)*rho_m11prime); 
end

% END OF BOUNDARY CONDITIONS - BOTTOM RIGHT WALL

% REGULARIZED BOUNDARY CONDITIONS - TOP LEFT WALL 
rhos = 0;
rhos_m20Istar = 0;
rhos_m11Istar = 0;
rhos_m02Istar = 0;

for i=1:9
   if (cx(i)<=0 && cy(i)>=0)
      rhos = rhos + fIn(i,1,ly);      
      rhos_m20Istar = rhos_m20Istar + fIn(i,1,ly)*(cx(i)*cx(i)-1/3); 
      rhos_m11Istar = rhos_m11Istar + fIn(i,1,ly)*(cx(i)*cy(i)); 
      rhos_m02Istar = rhos_m02Istar + fIn(i,1,ly)*(cy(i)*cy(i)-1/3);
   end
end

rhoprime = -12/(16+9*omega)*(-3*rhos-3*rhos_m02Istar+omega*3*rhos_m02Istar-...
            7*rhos_m11Istar+omega*7*rhos_m11Istar-3*rhos_m20Istar+...
            omega*3*rhos_m20Istar);
rho_m20prime = 2/9*(rhoprime+6*rhos_m11Istar+9*rhos_m20Istar);
rho_m11prime = 1/27*(7*rhoprime+18*rhos_m02Istar+132*rhos_m11Istar+18*rhos_m20Istar);
rho_m02prime = 2/9*(rhoprime+6*rhos_m11Istar+9*rhos_m02Istar);

for i=1:9 
    fIn(i,1,ly) =  w(i) * ... 
        ( rhoprime + ...
        9/2*(cx(i)*cx(i)-1/3)*rho_m20prime + ...
        9/2*(cy(i)*cy(i)-1/3)*rho_m02prime + ...    
        9*cx(i)*cy(i)*rho_m11prime);  
end
% END OF BOUNDARY CONDITIONS - TOP LEFT WALL

% REGULARIZED BOUNDARY CONDITIONS - TOP RIGHT WALL 
rhos = 0;
rhos_m20Istar = 0;
rhos_m11Istar = 0;
rhos_m02Istar = 0;

for i=1:9
   if (cx(i)>=0 && cy(i)>=0)
      rhos = rhos + fIn(i,lx,ly);      
      rhos_m20Istar = rhos_m20Istar + fIn(i,lx,ly)*(cx(i)*cx(i)-1/3); 
      rhos_m11Istar = rhos_m11Istar + fIn(i,lx,ly)*(cx(i)*cy(i)); 
      rhos_m02Istar = rhos_m02Istar + fIn(i,lx,ly)*(cy(i)*cy(i)-1/3);
   end
end

rhoprime = -12/(16+9*omega)*(-3*rhos-3*rhos_m02Istar+omega*3*rhos_m02Istar+...
            7*rhos_m11Istar-omega*7*rhos_m11Istar-3*rhos_m20Istar+...
            omega*3*rhos_m20Istar);
rho_m20prime = 2/9*(rhoprime-6*rhos_m11Istar+9*rhos_m20Istar);
rho_m11prime = 1/27*(-7*rhoprime-18*rhos_m02Istar+132*rhos_m11Istar-18*rhos_m20Istar);
rho_m02prime = 2/9*(rhoprime-6*rhos_m11Istar+9*rhos_m02Istar );
        

for i=1:9 
    fIn(i,lx,ly) =  w(i) * ... 
        ( rhoprime + ...
        9/2*(cx(i)*cx(i)-1/3)*rho_m20prime + ...
        9/2*(cy(i)*cy(i)-1/3)*rho_m02prime + ...    
        9*cx(i)*cy(i)*rho_m11prime); 
end

% END OF BOUNDARY CONDITIONS - TOP RIGHT WALL

% REGULARIZED BOUNDARY CONDITIONS - LEFT WALL 
rhos = zeros(1,1,ly-2);
rhos_m20Istar = zeros(1,1,ly-2);
rhos_m11Istar = zeros(1,1,ly-2);
rhos_m02Istar = zeros(1,1,ly-2);

for i=1:9
   if cx(i)<=0 
      rhos = rhos + fIn(i,1,2:(ly-1));      
      rhos_m20Istar = rhos_m20Istar + fIn(i,1,2:(ly-1))*(cx(i)*cx(i)-1/3); 
      rhos_m11Istar = rhos_m11Istar + fIn(i,1,2:(ly-1))*(cx(i)*cy(i)); 
      rhos_m02Istar = rhos_m02Istar + fIn(i,1,2:(ly-1))*(cy(i)*cy(i)-1/3);
   end
end

rhoprime = (12*rhos+(1-omega)*9*rhos_m20Istar)/(9+omega);
rho_m20prime = (2*rhos+15*rhos_m20Istar)/(9+omega);
rho_m11prime = 2*rhos_m11Istar;
rho_m02prime = 6/5*rhos_m02Istar;

for i=1:9 
    fIn(i,1,2:(ly-1)) =  w(i) * ... 
        ( rhoprime + ...
        9/2*(cx(i)*cx(i)-1/3)*rho_m20prime + ...
        9/2*(cy(i)*cy(i)-1/3)*rho_m02prime + ...    
        9*cx(i)*cy(i)*rho_m11prime); 

   end
% END OF BOUNDARY CONDITIONS - LEFT WALL

% REGULARIZED BOUNDARY CONDITIONS - RIGHT WALL 
rhos = zeros(1,1,ly-2);
rhos_m20Istar = zeros(1,1,ly-2);
rhos_m11Istar = zeros(1,1,ly-2);
rhos_m02Istar = zeros(1,1,ly-2);

for i=1:9
   if cx(i)>=0 
      rhos = rhos + fIn(i,lx,2:(ly-1));      
      rhos_m20Istar = rhos_m20Istar + fIn(i,lx,2:(ly-1))*(cx(i)*cx(i)-1/3); 
      rhos_m11Istar = rhos_m11Istar + fIn(i,lx,2:(ly-1))*(cx(i)*cy(i)); 
      rhos_m02Istar = rhos_m02Istar + fIn(i,lx,2:(ly-1))*(cy(i)*cy(i)-1/3);
   end
end

rhoprime = (12*rhos+(1-omega)*9*rhos_m20Istar)/(9+omega);
rho_m20prime = (2*rhos+15*rhos_m20Istar)/(9+omega);
rho_m11prime = 2*rhos_m11Istar;
rho_m02prime = 6/5*rhos_m02Istar;

for i=1:9 
    fIn(i,lx,2:(ly-1)) =  w(i) * ... 
        ( rhoprime + ...
        9/2*(cx(i)*cx(i)-1/3)*rho_m20prime + ...
        9/2*(cy(i)*cy(i)-1/3)*rho_m02prime + ...    
        9*cx(i)*cy(i)*rho_m11prime); 
end
% END OF BOUNDARY CONDITIONS - RIGHT WALL

 
rho = sum(fIn); 

% tau_xx = reshape ( ((cx.*cx-1/3)*reshape(fIn,9,lx*ly)), 1,lx,ly)-rho.*ux.*ux;
% tau_xy = reshape ( ((cx.*cy)*reshape(fIn,9,lx*ly)), 1,lx,ly)-rho.*ux.*uy;
% tau_yy = reshape ( ((cy.*cy-1/3)*reshape(fIn,9,lx*ly)), 1,lx,ly)-rho.*uy.*uy;

% COLLISION STEP 
 
for i=1:9 
    cu = 3*(cx(i)*ux+cy(i)*uy);     
    fEq(i,:,:) = rho .* w(i) .* ... 
        ( 1 + cu + 1/2*(cu.*cu) - ...
        3/2*(ux.^2+uy.^2) + ...
        27/2*ux.^2.*uy*(cx(i)^2-1/3)*cy(i)+ ...
        27/2*uy.^2.*ux*(cy(i)^2-1/3)*cx(i)+...
        81/4*(ux.^2.*uy.^2)*(cx(i)^2-1/3)*(cy(i)^2-1/3)); 
      fOut(i,:,:) = fIn(i,:,:)-omega*(fIn(i,:,:)-fEq(i,:,:));
    % fOut(i,:,:) = fEq(i,:,:)+(1-omega)*fnEq(i,:,:);  
end

tau_xx = reshape ( ((cx.*cx-1/3)*reshape(fOut,9,lx*ly)), 1,lx,ly)-rho.*ux.*ux;
tau_xy = reshape ( ((cx.*cy)*reshape(fOut,9,lx*ly)), 1,lx,ly)-rho.*ux.*uy;
tau_yy = reshape ( ((cy.*cy-1/3)*reshape(fOut,9,lx*ly)), 1,lx,ly)-rho.*uy.*uy;


% m20star = reshape ( ((cx.*cx -1/3)* reshape(fOut,9,lx*ly)), 1,lx,ly ) ./rho;
% m11star = reshape ( ((cx.*cy)* reshape(fOut,9,lx*ly)), 1,lx,ly ) ./rho;
% m02star = reshape ( ((cy.*cy -1/3)* reshape(fOut,9,lx*ly)), 1,lx,ly ) ./rho;
% 

% REGULARIZATION PROCEDURE
for i=1:9 
    cu = 3*(cx(i)*ux+cy(i)*uy);     

    fEq(i,:,:) = rho .* w(i) .* ... 
        ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2) + ...
        27/2*ux.^2.*uy*(cx(i)^2-1/3)*cy(i)+ ...
        27/2*uy.^2.*ux*(cy(i)^2-1/3)*cx(i)+...
        81/4*(ux.^2.*uy.^2)*(cx(i)^2-1/3)*(cy(i)^2-1/3)); 
    fnEq(i,:,:) = w(i) .* (...
        9/2 * (cx(i)^2 - 1/3) * tau_xx + ...
        9/2 * (cy(i)^2 - 1/3) * tau_yy + ...
        9 * (cx(i) * cy(i)) * tau_xy + ...
        27/2 * (uy .* tau_xx + 2 * ux .* tau_xy) .* (cx(i)^2 - 1/3) .* cy(i) + ...
        27/2 * (ux .* tau_yy + 2 * uy .* tau_xy) .* (cy(i)^2 - 1/3) .* cx(i) + ...
        81/4 * (uy.^2 .* tau_xx + 4 * ux .* uy .* tau_xy + uy.^2.*tau_yy) .* ...
        (cx(i)^2 - 1/3) .* (cy(i)^2 - 1/3));
    fOut(i,:,:) = fEq(i,:,:)+fnEq(i,:,:);
end


% STREAMING STEP 
for i=1:9 
    fIn(i,:,: ) = circshift(fOut(i,:,: ), [0,cx(i),cy(i)]); 
end

ec(cycle) = .5*sum(sum(ux.^2+uy.^2))/(lx*ly*uLid*uLid);

% % VISUALIZATION
 if (mod(cycle,tPlot)==0)
    cycle %#ok<NOPTS>

     ux_plot=reshape(ux,lx,ly);
     uy_plot=reshape(uy,lx,ly);
     grid on
     plot(ux_plot((lx+1)/2,:)/uLid/2,[0:ly-1]/(ly-1)-.5,'-.');
     hold on 
     plot([0:lx-1]/(lx-1)-.5,uy_plot(:,(ly+1)/2)/uLid/2,'-.');
     % % Plotar os dados de referência de Ghia
     plot(u_ghia / 2, y_ghia - .5, '-o');  % Velocidade u ao longo de y
     plot(x_ghia - .5, v_ghia / 2, '-o');  % Velocidade v ao longo de x
     % scatter(u_ghia / 2, y_ghia, 7, 'black', 'filled');  % Velocidade u 
     % scatter(x_ghia - .5, v_ghia, 7, 'red', 'filled');  % Velocidade v
     title('Velocidade nas Linhas Centrais')
     axis([-0.5 0.5 -0.5 0.5]);
     axis square      
     hold off

    % Configurações do gráfico
    title('Velocidade nas Linhas Centrais');
    xlabel('Coordenada Normalizada');
    ylabel('Velocidade Normalizada');
    drawnow     
 end
end

toc

plot(ec)

save ldcRe1e4L065hor
savefig('plotghia')

%u = reshape(sqrt(ux.^2+uy.^2),lx,ly); % para 
%imagesc(u(:,ly:-1:1)'./uLid);
%colorbar
%axis equal off; drawnow