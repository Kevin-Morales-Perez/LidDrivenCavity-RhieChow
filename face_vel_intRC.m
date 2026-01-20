function [u_f,v_f] = face_vel_intRC(u_f,v_f,u,v,A_P,p,dx,dy,nx,ny)
    
    %West - East faces ##################################################
    
    %West Boundary
    u_f(:,1)=0;

    %East faces of west boundary
    j=2;
    for i=1:ny
        u_f(i,j)=0.5*(u(i,j-1) + u(i,j)) + ...
            0.5*((1/A_P(i,j-1)) + (1/A_P(i,j)))*dy*(p(i,j-1) -p(i,j)) ...
            -0.25*dy*((1/A_P(i,j-1))*(p(i,j) -p(i,j-1))  + (1/A_P(i,j))*(p(i,j+1)  - p(i,j-1)));
            
    end

    %Interior Faces

    for i=1:ny
        for j=3:nx-1
            u_f(i,j)=0.5*(u(i,j-1) + u(i,j)) + ...
                0.5*((1/A_P(i,j-1)) + (1/A_P(i,j)))*dy*(p(i,j-1) -p(i,j))...
                -0.25*dy*((1/A_P(i,j-1))*(p(i,j) -p(i,j-2))  + (1/A_P(i,j))*(p(i,j+1)  - p(i,j-1)));
        end
    end

    j=nx;
    for i=1:ny
        u_f(i,j)=0.5*(u(i,j-1) + u(i,j)) + ...
            0.5*((1/A_P(i,j-1)) + (1/A_P(i,j)))*dy*(p(i,j-1) -p(i,j))...
            -0.25*dy*((1/A_P(i,j-1))*(p(i,j) -p(i,j-2))  + ...
            (1/A_P(i,j))*(p(i,j)  - p(i,j-1)));

    end

    %East Boundary      *************************************************
    u_f(:,nx+1)=0;%Solid wall


    %North - South faces %################################################
    
    %North Boundary     *************************************************
    v_f(1,:)=0;

    %South faces of north boundary
    i=2;
    for j=1:nx
        v_f(i,j)=0.5*(v(i-1,j) + v(i,j)) + ...
            0.5*((1/A_P(i,j)) + (1/A_P(i-1,j)))*dx*(p(i,j)-p(i-1,j))...
            -0.25*dx*((1/A_P(i,j))*(p(i-1,j)-p(i+1,j))  + ...
            (1/A_P(i-1,j))*(p(i-1,j)-p(i,j)));
    end


    %Interior Faces     *************************************************
    for i=3:ny-1
        for j=1:nx
            v_f(i,j)=0.5*(v(i-1,j) + v(i,j)) + ...
                0.5*((1/A_P(i,j)) + (1/A_P(i-1,j)))*dx*(p(i,j)-p(i-1,j))...
                -0.25*dx*((1/A_P(i,j))*(p(i-1,j)-p(i+1,j))  + ...
                (1/A_P(i-1,j))*(p(i-2,j)-p(i,j)));
        end
    end

    i=ny;
    for j=1:nx
        v_f(i,j)=0.5*(v(i-1,j) + v(i,j)) + ...
                0.5*((1/A_P(i,j)) + (1/A_P(i-1,j)))*dx*(p(i,j)-p(i-1,j))...
                -0.25*dx*((1/A_P(i,j))*(p(i-1,j)-p(i,j))  + ...
                (1/A_P(i-1,j))*(p(i-2,j)-p(i,j)));
    end

    %North faces of south boundary

    %South Boundary     **********************************************
    v_f(ny+1,nx)=0;

end

%************************************************************************