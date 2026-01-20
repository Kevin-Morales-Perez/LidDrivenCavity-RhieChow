function [A_W,A_N,A_E,A_S,A_P,Su_x,Su_y] = momentum_link_coeff(rho,mu,u0...
    ,dx,dy,u_f,v_f,p,nx,ny,A_W,A_N,A_E,A_S,A_P,Su_x,Su_y)
    
    %Momentum link coeffitiens

    Dh=(dy/dx)*mu;%diffusion flux for W-E nodes
    Dv=(dx/dy)*mu;%Difussion flux for S-N nodes

    %Central cells

    for i=2:ny-1
        for j=2:nx-1

            %compute mass fluxes
            fW=-rho*dy*u_f(i,j);
            fN=rho*dx*v_f(i,j);
            fE=rho*dy*u_f(i,j+1);
            fS=-rho*dx*v_f(i+1,j);

            massB = fW +  fN +  fE + fS;%Mass balance

            %Diagonal coeffitients
            A_W(i,j)=Dh + max(0,-fW);
            A_N(i,j)=Dv + max(0,-fN);
            A_E(i,j)=Dh + max(0,-fE);
            A_S(i,j)=Dv + max(0,-fS);

            %Central coeffitient
            A_P(i,j)= A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

            %Sources
            Su_x(i,j)=0.5*dy*(p(i,j-1) - p(i,j+1));
            Su_y(i,j)=0.5*dx*(p(i+1,j) - p(i-1,j));

        end
    end

    %EDGES
    %West
    j=1;
    for i=2:ny-1
        %compute mass fluxes
        %fW=-rho*dy*u_f(i,j);
        fN=rho*dx*v_f(i,j);
        fE=rho*dy*u_f(i,j+1);
        fS=-rho*dx*v_f(i+1,j);

        massB = fN +  fE + fS;%Mass balance

        %Diagonal coeffitients
        A_W(i,j)=2*Dh;
        A_N(i,j)=Dv + max(0,-fN);
        A_E(i,j)=Dh + max(0,-fE);
        A_S(i,j)=Dv + max(0,-fS);

        %Central coeffitient
        A_P(i,j)=A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

        %Sources
        Su_x(i,j)=0.5*dy*(p(i,j) - p(i,j+1));
        Su_y(i,j)=0.5*dx*(p(i-1,j) - p(i+1,j));

    end

    %North
    i=1;
    for j=2:nx-1

        %compute mass fluxes
        fW=-rho*dy*u_f(i,j);
        fE=rho*dy*u_f(i,j+1);
        fS=-rho*dx*v_f(i+1,j);

        massB = fW + fE + fS;%Mass balance

        %Diagonal coeffitients
        A_W(i,j)=Dh + max(0,-fW);
        A_N(i,j)=2*Dv;
        A_E(i,j)=Dh + max(0,-fE);
        A_S(i,j)=Dv + max(0,-fS);

        %Central coeffitient
        A_P(i,j)= A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

        %Sources
        Su_x(i,j)=0.5*dy*(p(i,j-1) - p(i,j+1)) + u0*A_N(i,j);
        Su_y(i,j)=0.5*dx*(p(i+1,j) - p(i,j));
        
    end

    %East
    j=nx;
    for i=2:ny-1

        %compute mass fluxes
        fW=-rho*dy*u_f(i,j);
        fN=rho*dx*v_f(i,j);
        fS=-rho*dx*v_f(i+1,j);

        massB = fW +  fN + fS;%Mass balance

        %Diagonal coeffitients
        A_W(i,j)=Dh + max(0,-fW);
        A_N(i,j)=Dv + max(0,-fN);
        A_E(i,j)=2*Dh;
        A_S(i,j)=Dv + max(0,-fS);

        %Central coeffitient
        A_P(i,j)= A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

        %Sources
        Su_x(i,j)=0.5*dy*(p(i,j-1) - p(i,j));
        Su_y(i,j)=0.5*dx*(p(i-1,j) - p(i+1,j));
    end

    %South
    i=ny;
    for j=2:nx-1

        %compute mass fluxes
        fW=-rho*dy*u_f(i,j);
        fN=rho*dx*v_f(i,j);
        fE=rho*dy*u_f(i,j+1);
        
        massB = fW +  fN +  fE ;%Mass balance

        %Diagonal coeffitients
        A_W(i,j)=Dh + max(0,-fW);
        A_N(i,j)=Dv + max(0,-fN);
        A_E(i,j)=Dh + max(0,-fE);
        A_S(i,j)=2*Dv;

        %Central coeffitient
        A_P(i,j)= A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

        %Sources
        Su_x(i,j)=0.5*dy*(p(i,j-1) - p(i,j+1));
        Su_y(i,j)=0.5*dx*(p(i,j) - p(i-1,j));

    end


    %Corners

    %West-North
    i=1;
    j=1;

    %compute mass fluxes
    fE=rho*dy*u_f(i,j+1);
    fS=-rho*dx*v_f(i+1,j);

    massB = fE + fS;%Mass balance

    %Diagonal coeffitients
    A_W(i,j)=2*Dh;
    A_N(i,j)=2*Dv;
    A_E(i,j)=Dh + max(0,-fE);
    A_S(i,j)=Dv + max(0,-fS);

    %Central coeffitient
    A_P(i,j)= A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

    %Sources
    Su_x(i,j)=0.5*dy*(p(i,j) - p(i,j+1)) + u0*A_N(i,j) ;
    Su_y(i,j)=0.5*dx*(p(i+1,j) - p(i,j));


    %East-North
    i=1;
    j=nx;

    %compute mass fluxes
    fW=-rho*dy*u_f(i,j);
    fS=-rho*dx*v_f(i+1,j);

    massB = fW +  fN +  fE + fS;%Mass balance

    %Diagonal coeffitients
    A_W(i,j)=Dh + max(0,-fW);
    A_N(i,j)=Dv;
    A_E(i,j)=Dh;
    A_S(i,j)=Dv + max(0,-fS);

    %Central coeffitient
    A_P(i,j)= A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

    %Sources

    Su_x(i,j)=0.5*dy*(p(i,j-1) - p(i,j)) + u0*A_N(i,j);
    Su_y(i,j)=0.5*dx*(p(i+1,j) - p(i,j));

    %East-South
    i=ny;
    j=nx;
    
    %compute mass fluxes
    fW=-rho*dy*u_f(i,j);
    fN=rho*dx*v_f(i,j);

    massB = fW +  fN ;%Mass balance

    %Diagonal coeffitients
    A_W(i,j)=Dh + max(0,-fW);
    A_N(i,j)=Dv + max(0,-fN);
    A_E(i,j)=2*Dh;
    A_S(i,j)=2*Dv;

    %Central coeffitient
    A_P(i,j)= A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

    %Sources
    Su_x(i,j)=0.5*dy*(p(i,j-1) - p(i,j));
    Su_y(i,j)=0.5*dx*(p(i,j) - p(i-1,j));

    %West-South
    i=ny;
    j=1;

    %compute mass fluxes
    fN=rho*dx*v_f(i,j);
    fE=rho*dy*u_f(i,j+1);

    massB = fN +  fE;%Mass balance

    %Diagonal coeffitients
    A_W(i,j)=2*Dh;
    A_N(i,j)=Dv + max(0,-fN);
    A_E(i,j)=Dh + max(0,-fE);
    A_S(i,j)=2*Dv;

    %Central coeffitient
    A_P(i,j)= A_W(i,j) + A_N(i,j) + A_E(i,j) + A_S(i,j) + massB;

    %Sources
    Su_x(i,j)=0.5*dy*(p(i,j) - p(i,j+1));
    Su_y(i,j)=0.5*dx*(p(i,j) - p(i-1,j));

end

