% Lid driven cavity in co-located mesh

%WRITTEN BY: KEVIN MORALES 
%AERONAUTICAL ENGINEERING - INSTITUTO POLITECNICO NACIONAL - CIUDAD DE
%MEXICO
%-kevin27182mora@gmail.com
%% DATA

%physical parameters
%Water at 20°
rho=998; %density
mu=1.003e-3;%kinematic viscosity
nu=mu/rho;%Dynamic viscosity
u0=10*0.0100501;% xvelocity of the lid 

%Domain
lenght_cvty=0.1;
height_cvty=0.1;

%Reynolds Number
Re=rho*u0*(0.5*(lenght_cvty  + height_cvty ))/mu;
disp(Re)
%% MESH
nx=70;
ny=70;
ncells=nx*ny;

dx=lenght_cvty/nx;%x grid spacing 
dy=height_cvty/ny;%y grid spacing 

%Nodal Point positions 
x_centers=0.5*dx:dx:lenght_cvty-0.5*dx;
y_centers=height_cvty-0.5*dy:-dy:0.5*dy;

[Xctrs,Yctrs]=meshgrid(x_centers,y_centers);

%plot(X,Y,'o')% check the positions of each point in the grid 
%axis equal


%% VARIABLES
u=zeros(ny,nx); %velocity in X axis
v=zeros(ny,nx); %Velocity in Y axis
p=ones(ny,nx); %Pressure
%p_prime=zeros(ny,nx);% Pressure correction

u_f=zeros(ny,nx+1); %x velocity at faces
v_f=zeros(ny+1,nx); %y velocity at faces

u_f(:,1)=u0; %Imposing Inlet velocity

%Coeffitients for momentum equations 
%Neighborhood coeffitients 
A_W=zeros(ny,nx);
A_N=zeros(ny,nx);
A_E=zeros(ny,nx);
A_S=zeros(ny,nx);

%Central coeffitient
A_P=zeros(ny,nx);

%Sources 
Su_x=zeros(ny,nx);
Su_y=zeros(ny,nx);

%Coeffitients for pressure correction equations 
%Neighborhood coeffitients 
Ap_W=zeros(ny,nx);
Ap_N=zeros(ny,nx);
Ap_E=zeros(ny,nx);
Ap_S=zeros(ny,nx);

%Central coeffitient
Ap_P=zeros(ny,nx);

%Sources 
Su_p=zeros(ny,nx);

%Underelaxation factors
alpha_uv=0.7;  %<---------------------     x-y momentum 
alpha_p=0.005; %<---------------------     pressure correction

%Tolerance for inner iterations 
epsilon_u=1e-21;
epsilon_v=1e-21;
epsilon_p= 1e-21;

%Raw Residual
rsid_x=zeros(ny,nx); %X momentum
rsid_y=zeros(ny,nx); %Y momentum
rsid_p=zeros(ny,nx);% Pressure correction eq.
rsid_cont=zeros(ny,nx); % Continuity

%error from residual
err_x=1;
err_y=1;
err_p=1;
err_cont=1;

max_iterations=3000;% Max outer iterations <---------------------
max_iterations_u=5;% Max iterations for momentum eqs.
max_iterations_v=5;% Max iterations for momentum  y eq .
max_iterations_p=40;% Max iterations for pressure eq.
residualsMat=zeros(max_iterations,4);% Residual Matrix
inneritcontMat=zeros(max_iterations,3);% Inner iterations counter matrix 
gmrsErrors=zeros(max_iterations,3);%Comment if GMRES not used, add more 
% columns if GMRES is used in more than 1 solver
error_tgt=1e-20; %Target error 
max_residual=1e10; % unstable value flag
convergedFlg=false; %Flag for convergence

iterations_cont=0;

%%  INDEX CORRESPONDENCE BETWEEN MATRIX FORM AND VECTOR FORM

indx_mat=zeros(ncells,2);% MAtrix indexes
indx_k=zeros(ncells,1);% Vector indexes 

%k=(i-1)*nx + j

%P(center)  (i,j)  --->  k
%West       (i,j-1) ---> k - 1
%North      (i-1,j) ---> k - nx
%East       (i,j+1) ---> k + 1 
%South      (i+1,j) ---> k + 1

for i=1:ny
    for j=1:nx
        k=(i-1)*nx + j;
        indx_k(k)=k;
        indx_mat(k,:)=[i,j];
    end
end

%Interior, Edges , Corner Indexes

%vector indexes
indx_cor_WN=1;%West North corner (1)
indx_edg_N=2:nx-1;%North edge (2)
indx_cor_EN=nx;%East North corner(3)
indx_edg_W=nx+1:nx:(ny-2)*nx+1;%West Edge (4)
indx_edg_E=nx*2:nx:(ny-1)*nx;%East Edge (5)
indx_cor_WS=(ny-1)*nx +1;%West South corner (6)
indx_edg_S=(ny-1)*nx +2:ncells-1;%South Edge (7)
indx_cor_ES=ncells;%East South corner (8)
indx_interior=zeros(ny-2,nx-2);%indexes for interior cells

for i=2:ny-1
    indx_interior(i-1,:)=(i-1)*nx+2:1:i*nx-1;
end
indx_interior=reshape(transpose(indx_interior),1,[]);


%% GENERALIZED  MINIMUM  RESIDUAL METHOD SETTINGS
%x momentum
A_xsp=sparse(ncells,ncells);%ALL coeffitients
suX_vec=zeros(ncells,1);%Source Vector
restart_x=1;%Restart vector space
tol_x=1e-3;%tolerance desired
maxit_x=5;%Inner it

err_x2=1;%Output error from GRMES

%y momentum
A_ysp=sparse(ncells,ncells);%ALL coeffitients
suY_vec=zeros(ncells,1);%Source Vector
restart_y=1;%Restart vector space
tol_y=1e-3;%tolerance desired
maxit_y=5;%Inner it

err_y2=1;%Output error from GRMES


%pressure correction
A_psp=sparse(ncells,ncells);%ALL coeffitients
suP_vec=zeros(ncells,1);%Source Vector
restart_p=2;%Restart vector space
tol_p=1e-15;%tolerance desired
maxit_p=4;%Inner it

err_p2=1;%Output error from GRMES


%% SOLVER SETTINGS
%G-S:GAUSS SEIDEL
%GMRES:GENERALIZED MINIMUM RESIDUAL METHOD
solver_x_momentum="G-S";
solver_y_momentum="G-S";
solver_p_correction="GMRES";

%% SOLVER
tic
while convergedFlg==false
    
    iterations_cont=iterations_cont+1;

    % #1 MOMENTUM LINK COEFFITIENTS
    [A_W,A_N,A_E,A_S,A_P,Su_x,Su_y] = momentum_link_coeff(rho,mu,u0...
    ,dx,dy,u_f,v_f,p,nx,ny,A_W,A_N,A_E,A_S,A_P,Su_x,Su_y);

    % #2 SOLVE  X MOMENTUM

    switch solver_x_momentum
        case "G-S"
    
            [u,rsid_x,err_x,it_innx] = GaussSeidel_GeneralSolver(A_P,A_W,...
            A_N,A_E,A_S,Su_x,rsid_x,epsilon_u,err_x,max_iterations_u,...
            alpha_uv,nx,ny,u,u);
        case "GMRES"
            [u,rsid_x,err_x,err_x2,int_innx] = gmresSolver(A_P,A_W,A_N,A_E,A_S,Su_x,...
            A_xsp,suX_vec,restart_x,tol_x,maxit_x,...
            indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,indx_edg_S,...
            indx_cor_WN,indx_cor_EN,indx_cor_ES,indx_cor_WS, ...
            indx_mat,ny,nx,ncells,u,...
            alpha_uv);
    end
     
    %#3 SOLVE  Y MOMENTUM

    switch solver_y_momentum
        case "G-S"
       
            [v,rsid_y,err_y,it_inny] = GaussSeidel_GeneralSolver(A_P,A_W,...
            A_N,A_E,A_S,Su_y,rsid_y,epsilon_v,err_y,max_iterations_v,...
            alpha_uv,nx,ny,v,v);

        case "GMRES"
            [v,rsid_y,err_y,err_y2,int_inny] = gmresSolver(A_P,A_W,A_N,A_E,A_S,Su_y,...
            A_ysp,suY_vec,restart_y,tol_y,maxit_y,...
            indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,indx_edg_S,...
            indx_cor_WN,indx_cor_EN,indx_cor_ES,indx_cor_WS, ...
            indx_mat,ny,nx,ncells,v,...
            alpha_uv);

    end

    %Apply underelaxation to main diagonal coeffitients
    A_P=A_P/alpha_uv;

    % #4 FACE VELOCITY COMPUTATION USING RHIE-CHOW INTERPOLATION
    [u_f,v_f] = face_vel_intRC(u_f,v_f,u,v,A_P,p,dx,dy,nx,ny);
 
    % #5 PRESSURE CORRECTION LINK COEFFITIENTS
    [Ap_W,Ap_N,Ap_E,Ap_S,Ap_P,Su_p] = pressCorr_link_coeff(u_f,v_f,...
    dx,dy,A_P,Ap_W,Ap_N,Ap_E,Ap_S,Ap_P,Su_p,nx,ny);

    % #6 SOLVE PRESSURE CORRECTION

    %restart pressure correction
    p_prime=zeros(ny,nx);
     switch solver_p_correction
        case "G-S"
        
            [p_prime,rsid_p,err_p,int_innp]=GaussSeidel_GeneralSolver(Ap_P,Ap_W,...
            Ap_N,Ap_E,Ap_S,Su_p,rsid_p,epsilon_p,err_p,max_iterations_p,...
            1,nx,ny,p_prime,p_prime);

        case "GMRES"
            [p_prime,rsid_p,err_p,err_p2,int_innp] = gmresSolver(Ap_P,Ap_W,Ap_N,Ap_E,Ap_S,Su_p,...
            A_psp,suP_vec,restart_p,tol_p,maxit_p,...
            indx_interior,indx_edg_W,indx_edg_N,indx_edg_E,indx_edg_S,...
            indx_cor_WN,indx_cor_EN,indx_cor_ES,indx_cor_WS, ...
            indx_mat,ny,nx,ncells,p_prime,...
            1);

     end

    % #7 CORRECT PRESSURE

    p =p + alpha_p*p_prime;

    % #8 CORRECT FACE VELOCITY
    [u_f,v_f] = fvel_correct(u_f,v_f,p_prime,A_P,1,dx,dy,nx,ny);

    % #9 CORRECT CELL CENTER VELOCITY
    [u,v] = cvel_correct(u,v,p_prime,A_P,1,dx,dy,nx,ny);

    % #10 CHECK CONVERGENCE CRITERIA 
    %Check continuity
    [rsid_cont,err_cont] = checkContinuity(u_f,v_f,rsid_cont,dx,dy,nx,ny);

    residualsMat(iterations_cont,:)=[err_x,err_y,err_p,err_cont];
    inneritcontMat(iterations_cont,:)=[it_innx,it_inny,int_innp];
    gmrsErrors(iterations_cont,:)=[err_x2,err_y2,err_p2];
    disp(residualsMat(iterations_cont,:)) %display current error


    if err_x < error_tgt && err_y < error_tgt && err_cont < error_tgt ...
            && err_p < error_tgt
       convergedFlg=true;
       fprintf("Converged  at iteration\n")
       disp(iterations_cont)
       
    elseif err_x>max_residual || err_y > max_residual ||...
            err_cont > max_residual || err_p>max_residual
        fprintf("Unstable solution iterations stopped at iteration ")
        disp(iterations_cont)
        break
    elseif iterations_cont >= max_iterations
        fprintf("Max iterations reached \n")
        disp(iterations_cont)
        break
    elseif isnan(err_x) || isnan(err_y) || isnan(err_cont) || isnan(err_p)
        fprintf("The system became undetermined  at iteration \n")
        disp(iterations_cont)
        break
    else
        convergedFlg=false;
    end

end
toc
%% POST PROCESS

%COMPUTATION OF OTHER VARIABLES OF INTEREST

%VELOCITY MAGNITUDE
vel= sqrt(u.^2 + v.^2); 


%RESIDUALS 
figure(1)
plot(1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,1),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,2),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,3),...
    1:iterations_cont-1,...
    residualsMat(1:iterations_cont-1,4))
legend("u vel","v vel","Pressure correction","Continuity")
title("Residuals")
xlabel("Iterations")
ylabel("Residual")
yscale log

figure(2)
plot(1:iterations_cont-1,...
    inneritcontMat(1:iterations_cont-1,1),...
    1:iterations_cont-1,...
    inneritcontMat(1:iterations_cont-1,2),...
    1:iterations_cont-1,...
    inneritcontMat(1:iterations_cont-1,3))
legend("u vel","v vel","Pressure correction")
title("Inner iterations per outer iteration")
xlabel("Iterations")
ylabel("Inner iterations")


figure(3)
contourf(Xctrs,Yctrs,p, 20, 'LineColor', 'none')
title("Pressure Distribution (N/m^2)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(4)
contourf(Xctrs,Yctrs,rsid_cont, 20, 'LineColor', 'none')
title("Continuity Residual (Kg*m^2/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal


figure(5)
contourf(Xctrs,Yctrs,u, 20, 'LineColor', 'none')
title("Velocity in x axis (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(6)
contourf(Xctrs,Yctrs,v, 20, 'LineColor', 'none')
title("Velocity in y axis (m/s)")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis equal

figure(7)
streamslice(Xctrs, Yctrs, u, v);
title("streamlines")
xlabel("Lenght (m)")
ylabel("Height (m)")
axis equal

figure(8)
contourf(Xctrs, Yctrs, vel, 20, 'LineColor', 'none')
colorbar
hold on 
quiver(Xctrs, Yctrs,u,v,'k')
title("Velocity vector field with magnitude m/s")
xlabel("Height (m)")
ylabel("Lenght (m)")
colormap jet
axis equal
hold off


