function [u_f,v_f] = fvel_correct(u_f_star,v_f_star,p_prime,A_P,alpha_uv,dx,dy,nx,ny)
%FACE VELOCITY CORRECTION WITH PRESSURE CORRECTION FIELD
    
    u_f=u_f_star;
    v_f=v_f_star;
    
    %West - East Faces
    for i=1:ny
        for j=2:nx
            u_f(i,j) = u_f_star(i,j) + alpha_uv*0.5*((1/A_P(i,j-1)) + ...
                (1/A_P(i,j)))*dy*(p_prime(i,j-1) -p_prime(i,j));
        end
    end

    %North - South Faces
    for i=2:ny
        for j=1:nx
            v_f(i,j)= v_f_star(i,j) + alpha_uv*0.5*((1/A_P(i-1,j)) + ...
                (1/A_P(i,j)))*dx*(p_prime(i,j) -p_prime(i-1,j)); 
        end
    end

end