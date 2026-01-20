function [u,v] = cvel_correct(u_star,v_star,p_prime,A_P,alpha_uv,dx,dy,nx,ny)
%CELL CENTER VELOCITY CORRECTION WITH PRESSURE CORRECTION FIELD
    u=u_star;
    v=v_star;
    
    %x axis
    %West Edge
    j=1;
    for i=1:ny
        u(i,j)=u_star(i,j) + alpha_uv*0.5*(dy/A_P(i,j))*(p_prime(i,j) - p_prime(i,j+1));
    end

    
    %nterior cells 
    for i=1:ny
        for j=2:nx-1
            u(i,j)=u_star(i,j) + alpha_uv*0.5*(dy/A_P(i,j))*(p_prime(i,j-1) - p_prime(i,j+1));
        end
    end

    %East Edge
    j=nx;
    for i=1:ny
        u(i,j)=u_star(i,j) + alpha_uv*0.5*(dy/A_P(i,j))*(p_prime(i,j-1) - p_prime(i,j));
    end

    %y axis
    %North Edge
    i=1;
    for j=1:nx
        v(i,j)=v_star(i,j) + alpha_uv*0.5*(dx/A_P(i,j))*(p_prime(i+1,j)-p_prime(i,j));
    end

    %Interior cells
    for i=2:ny-1
        for j=1:nx
            v(i,j)=v_star(i,j) + alpha_uv*0.5*(dx/A_P(i,j))*(p_prime(i+1,j)-p_prime(i-1,j));
        end
    end

    %South Edge
    i=ny;
    for j=1:nx
        v(i,j)=v_star(i,j) + alpha_uv*0.5*(dx/A_P(i,j))*(p_prime(i,j)-p_prime(i-1,j));
    end



end