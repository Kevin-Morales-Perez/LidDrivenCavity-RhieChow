function [Ap_W,Ap_N,Ap_E,Ap_S,Ap_P,Su_p] = pressCorr_link_coeff(u_f,v_f,...
    dx,dy,A_P,Ap_W,Ap_N,Ap_E,Ap_S,Ap_P,Su_p,nx,ny)
    
    %Pressure correction link coeffitients

    %Interior cells

    for i=2:ny-1
        for j=2:nx-1
            
            %Neighborhod coeffitiens 
            Ap_W(i,j)=0.5*(1/A_P(i,j-1) + 1/A_P(i,j))*(dy^2);
            Ap_N(i,j)=0.5*(1/A_P(i-1,j) + 1/A_P(i,j))*(dx^2);
            Ap_E(i,j)=0.5*(1/A_P(i,j+1) + 1/A_P(i,j))*(dy^2);
            Ap_S(i,j)=0.5*(1/A_P(i+1,j) + 1/A_P(i,j))*(dx^2);
            
            %Main Diagonal
            Ap_P(i,j)=Ap_W(i,j) + Ap_N(i,j) + Ap_E(i,j) + Ap_S(i,j);
            
            %Sources (Mass inbalance)
            Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
                dy*(v_f(i+1,j) - v_f(i,j));
                        
        end
    end

    %EDGES
    %West
    j=1;
    for i=2:ny-1

        %Neighborhod coeffitiens 
        Ap_N(i,j)=0.5*(1/A_P(i-1,j) + 1/A_P(i,j))*(dx^2);
        Ap_E(i,j)=0.5*(1/A_P(i,j+1) + 1/A_P(i,j))*(dy^2);
        Ap_S(i,j)=0.5*(1/A_P(i+1,j) + 1/A_P(i,j))*(dx^2);
        
        %Main Diagonal
        Ap_P(i,j)=Ap_N(i,j) + Ap_E(i,j) + Ap_S(i,j);
        
        %Sources (Mass inbalance)
        Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
            dy*(v_f(i+1,j) - v_f(i,j));


    end

    %North
    i=1;
    for j=2:nx-1
        %Neighborhod coeffitiens 
        Ap_W(i,j)=0.5*(1/A_P(i,j-1) + 1/A_P(i,j))*(dy^2);
        Ap_E(i,j)=0.5*(1/A_P(i,j+1) + 1/A_P(i,j))*(dy^2);
        Ap_S(i,j)=0.5*(1/A_P(i+1,j) + 1/A_P(i,j))*(dx^2);
        
        %Main Diagonal
        Ap_P(i,j)=Ap_W(i,j) + Ap_E(i,j) + Ap_S(i,j);
        
        %Sources (Mass inbalance)
        Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
            dy*(v_f(i+1,j) - v_f(i,j));
    end

    %East
    j=nx;
    for i=2:ny-1

        %Neighborhod coeffitiens 
        Ap_W(i,j)=0.5*(1/A_P(i,j-1) + 1/A_P(i,j))*(dy^2);
        Ap_N(i,j)=0.5*(1/A_P(i-1,j) + 1/A_P(i,j))*(dx^2);
        Ap_S(i,j)=0.5*(1/A_P(i+1,j) + 1/A_P(i,j))*(dx^2);
        
        %Main Diagonal
        Ap_P(i,j)=Ap_W(i,j) + Ap_N(i,j) + Ap_S(i,j);
        
        %Sources (Mass inbalance)
        Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
            dy*(v_f(i+1,j) - v_f(i,j));

    end

    %South
    i=ny;
    for j=2:nx-1

        %Neighborhod coeffitiens 
        Ap_W(i,j)=0.5*(1/A_P(i,j-1) + 1/A_P(i,j))*(dy^2);
        Ap_N(i,j)=0.5*(1/A_P(i-1,j) + 1/A_P(i,j))*(dx^2);
        Ap_E(i,j)=0.5*(1/A_P(i,j+1) + 1/A_P(i,j))*(dy^2);
        
        %Main Diagonal
        Ap_P(i,j)=Ap_W(i,j) + Ap_N(i,j) + Ap_E(i,j);
        
        %Sources (Mass inbalance)
        Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
            dy*(v_f(i+1,j) - v_f(i,j));

    end

 
    %CORNERS

    %West-North
    i=1;
    j=1;

    %Neighborhod coeffitiens 
    Ap_E(i,j)=0.5*(1/A_P(i,j+1) + 1/A_P(i,j))*(dy^2);
    Ap_S(i,j)=0.5*(1/A_P(i+1,j) + 1/A_P(i,j))*(dx^2);
    
    %Main Diagonal
    Ap_P(i,j)=Ap_E(i,j) + Ap_S(i,j);
    
    %Sources (Mass inbalance)
    Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
        dy*(v_f(i+1,j) - v_f(i,j));

    %East-North
    i=1;
    j=nx;

    %Neighborhod coeffitiens 
    Ap_W(i,j)=0.5*(1/A_P(i,j-1) + 1/A_P(i,j))*(dy^2);
    Ap_S(i,j)=0.5*(1/A_P(i+1,j) + 1/A_P(i,j))*(dx^2);
    
    %Main Diagonal
    Ap_P(i,j)=Ap_W(i,j) +  Ap_S(i,j);
    
    %Sources (Mass inbalance)
    Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
        dy*(v_f(i+1,j) - v_f(i,j));

    %East-South
    i=ny;
    j=nx;

    %Neighborhod coeffitiens 
    Ap_W(i,j)=0.5*(1/A_P(i,j-1) + 1/A_P(i,j))*(dy^2);
    Ap_N(i,j)=0.5*(1/A_P(i-1,j) + 1/A_P(i,j))*(dx^2);
    
    %Main Diagonal
    Ap_P(i,j)=Ap_W(i,j) + Ap_N(i,j);
    
    %Sources (Mass inbalance)
    Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
        dy*(v_f(i+1,j) - v_f(i,j));

    %West-South
    i=ny;
    j=1;

    %Neighborhod coeffitiens 
    Ap_N(i,j)=0.5*(1/A_P(i-1,j) + 1/A_P(i,j))*(dx^2);
    Ap_E(i,j)=0.5*(1/A_P(i,j+1) + 1/A_P(i,j))*(dy^2);
    
    %Main Diagonal
    Ap_P(i,j)=Ap_N(i,j) + Ap_E(i,j);
    
    %Sources (Mass inbalance)
    Su_p(i,j)=dx*(u_f(i,j) - u_f(i,j+1)) + ...
        dy*(v_f(i+1,j) - v_f(i,j));


end