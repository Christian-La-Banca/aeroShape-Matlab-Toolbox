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
        for i = 1:shape.surfacesArray(surf).nDivisions_i
            for j = 1:shape.surfacesArray(surf).nDivisions_j
                cosPhi(j,i) = shape.surfacesArray(surf).dS_vec{j,i}*Vvec_normalized/norm(shape.surfacesArray(surf).dS_vec{j,i});
    
                % Masking all values behind the body
                if cosPhi(j,i) > 0
                    cosPhi(j,i) = 0;
                end
            end
        end
        CP_currentSurface = options.Cpt2 * cosPhi.^2;
        shape.surfacesArray(surf).setCP(CP_currentSurface);

    end
    disp(['CP calculated using Newtoninan Hipersonic Flow in: ',num2str(toc*1000), ' miliseconds'])
end