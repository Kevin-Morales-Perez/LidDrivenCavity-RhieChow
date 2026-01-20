function [rsid_cont,err_cont] = checkContinuity(u_f,v_f,rsid_cont,...
    dx,dy,nx,ny)
    %CONTINUITY CHECK

    for i=1:ny
        for j=1:nx
            rsid_cont(i,j)=dx*(u_f(i,j+1) -u_f(i,j)) + dy*(v_f(i,j) - v_f(i+1,j));
        end
    end

    err_cont=rms(rsid_cont(:));


end