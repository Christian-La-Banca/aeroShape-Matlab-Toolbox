classdef shapeHipersonic < genericShape
    methods
        function solveCP(obj,alpha,beta)
            arguments
                obj shapeHipersonic
                alpha {mustBeReal}
                beta {mustBeReal}
            end

            Cpt2 = 2;

            Vvec_normalized = [cos(alpha)*cos(beta);cos(alpha)*sin(beta);sin(alpha)];

            for surf = 1:obj.numberOfSurfaces
                cosPhi = zeros(size(obj.surfacesArray(surf).X_data));
                for i = 1:obj.surfacesArray(surf).nDivisions_i
                    for j = 1:obj.surfacesArray(surf).nDivisions_j
                        cosPhi(j,i) = obj.surfacesArray(surf).dS_vec{j,i}*Vvec_normalized/norm(obj.surfacesArray(surf).dS_vec{j,i});

                        % Masking all values behind the body
                        if cosPhi(j,i) > 0
                            cosPhi(j,i) = 0;
                        end
                    end
                end
                CP_currentSurface = Cpt2 * cosPhi.^2;
                obj.surfacesArray(surf).setCP(CP_currentSurface);
            end
        end
    end
end