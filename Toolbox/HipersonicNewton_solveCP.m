function HipersonicNewton_solveCP(shape,alpha,beta,options)
arguments
    shape aeroShape
    alpha (1,1) {mustBeReal}
    beta (1,1) {mustBeReal}
    options.Cpt2 (1,1) {mustBeReal} = 2;
end

    Vvec_normalized = [cos(alpha)*cos(beta);cos(alpha)*sin(beta);sin(alpha)];
    tic
    for surf = 1:shape.numberOfSurfaces
        cosPhi = zeros(size(shape.surfacesArray(surf).X_data));
        dS_vec = {shape.surfacesArray(surf).dSx;shape.surfacesArray(surf).dSy;shape.surfacesArray(surf).dSz};
        for i = 1:shape.surfacesArray(surf).nDivisions_i
            for j = 1:shape.surfacesArray(surf).nDivisions_j
                dS = [dS_vec{1}(j,i),dS_vec{2}(j,i),dS_vec{3}(j,i)];
                cosPhi(j,i) = dS*Vvec_normalized/norm(dS);
    
                % Masking all values behind the body
                if cosPhi(j,i) > 0
                    cosPhi(j,i) = 0;
                end
            end
        end
        CP_currentSurface = options.Cpt2 * cosPhi.^2;
        CP_currentSurface(isnan(CP_currentSurface)) = 0;

        shape.surfacesArray(surf).setCP(CP_currentSurface);
    end
    disp(['CP calculated using Newtoninan Hipersonic Flow in: ',num2str(toc*1000), ' miliseconds'])
end