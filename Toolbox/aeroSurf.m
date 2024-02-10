classdef aeroSurf < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name;
    end

    properties (GetAccess = public,SetAccess = protected)
        X_data;
        Y_data;
        Z_data;
        dS_vec;
        Cp_data;
        nDivisions_i;
        nDivisions_j;
    end


    methods
        function obj = aeroSurf()
        end

        function setGeometry(obj,X,Y,Z,nameStr)
            arguments
                obj aeroSurf
                X {mustBeNumeric,mustBeReal,mustBeNonNan}
                Y {mustBeNumeric,mustBeReal,mustBeNonNan}
                Z {mustBeNumeric,mustBeReal,mustBeNonNan}
                nameStr {mustBeTextScalar} = 'aeroSurfaceDefault';
            end

            if all(size(X) ~= size(Y)) || all(size(X) ~= size(Z)) || all(size(Y) ~= size(Z))
                error('All 3 matrices for geometry data must be the same size')
            end

            obj.X_data = X;
            obj.Y_data = Y;
            obj.Z_data = Z;
            obj.name = nameStr;
            obj.nDivisions_i = size(X,2);
            obj.nDivisions_j = size(X,1);
            
            %Calculate normal vectors using cross product
            calculateNormals(obj);

            %Extrapolate 0 length normals
            extrapolateNullLengthNormals(obj);

            disp(['Normals calculated for surface: ', obj.name])
        end

        function flipNormals(obj)
            for i = 1:obj.nDivisions_i
                for j = 1:obj.nDivisions_j
                    obj.dS_vec{j,i} = -obj.dS_vec{j,i};
                end
            end
            disp(['Normals flipped for surface: ', obj.name])
        end

        function setCP(obj,CPmatrix)
            if all(size(CPmatrix) ~= size(obj.X_data))
                error('CP matrix must be the same size as geometry matrix')
            end

            obj.Cp_data = CPmatrix;
        end

        function [CF] = getForces(obj)

            if ~all(size(obj.Cp_data))
                error(['CP_data not populated for surface: ',obj.name,'. Use obj.setCP() to populate CP field'])
            end
            %Calculate integral as Riemmann sum
            CF = [0,0,0];
            for i = 1:obj.nDivisions_i
                for j = 1:obj.nDivisions_j
                    dS_face = obj.getFace_dS(j,i);
                    CP_face = obj.getFace_CP(j,i);
                    CF = CF - CP_face*dS_face;
                end
            end

            CF = transpose(CF);
        end

        function [CF,CM] = getForcesAndMoments(obj,momentsPoint)
            arguments
                obj aeroSurf
                momentsPoint (1,3)
            end

            if ~all(size(obj.Cp_data))
                error(['CP_data not populated for surface: ',obj.name,'. Use obj.setCP() to populate CP field'])
            end

            CF = [0,0,0];
            CM = [0,0,0];
            for i = 1:obj.nDivisions_i
                for j = 1:obj.nDivisions_j
                    dS_face = obj.getFace_dS(j,i);
                    CP_face = obj.getFace_CP(j,i);

                    CF_face = - CP_face*dS_face;
                    CF = CF + CF_face;

                    distanceVector = [obj.X_data(j,i),obj.Y_data(j,i),obj.Z_data(j,i)] - momentsPoint;
                    CM = CM + cross(distanceVector,CF_face);
                end
            end

            CF = transpose(CF);
            CM = transpose(CM);
        end

        function drawGeometry(obj)
            surf(obj.X_data,obj.Y_data,obj.Z_data,"EdgeColor",[128,128,128]./255,'EdgeAlpha',0.1)
            hold on
            axis equal
            zlabel('\bf Z')
            ylabel('\bf Y')
            xlabel('\bf X')
            colormap("turbo")            
        end

        function drawNormals(obj)
            drawGeometry(obj)

            for i = 1:(floor(obj.nDivisions_i/6)):obj.nDivisions_i
                for j = 1:(floor(obj.nDivisions_j/6)):obj.nDivisions_j
                    x = obj.X_data(j,i);
                    y = obj.Y_data(j,i);
                    z = obj.Z_data(j,i);
                    dS = 0.25*obj.dS_vec{j,i}/norm(obj.dS_vec{j,i});
                    arrow3([x,y,z],[x+dS(1),y+dS(2),z+dS(3)],'r-1',0.5,1)
                end
            end
        end

        function drawCP(obj)

            if ~all(size(obj.Cp_data))
                error(['CP_data not populated for surface: ',obj.name,'. Use obj.setCP() to populate CP field'])
            end

            surf(obj.X_data,obj.Y_data,obj.Z_data,obj.Cp_data,"EdgeColor",[128,128,128]./255,'EdgeAlpha',0.1)
            hold on
            axis equal
            zlabel('\bf Z')
            ylabel('\bf Y')
            xlabel('\bf X')
            colormap("turbo")
            cb = colorbar;
            ylabel(cb,'CP')
        end
    end

    methods (Access = protected)
        function calculateNormals(obj)
            obj.dS_vec = cell(size(obj.X_data,1),size(obj.X_data,2));

            verteces = {[0,0,0],[0,0,0],[0,0,0],[0,0,0]};
            for i = 1:obj.nDivisions_i
                for j=1:obj.nDivisions_j

                    if  i < obj.nDivisions_i && j < obj.nDivisions_j
                        queryIndeces = {[i,j],[i+1,j],[i,j+1],[i+1,j+1]};
                        for k = 1:4
                            currentVertexIndeces = queryIndeces{k};
                            verteces{k} = [obj.X_data(currentVertexIndeces(2),currentVertexIndeces(1)),...
                                           obj.Y_data(currentVertexIndeces(2),currentVertexIndeces(1)),...
                                           obj.Z_data(currentVertexIndeces(2),currentVertexIndeces(1))];
                        end

                        directionVector_i = verteces{2} - verteces{1};
                        directionVector_j = verteces{3} - verteces{1};

                    elseif i == obj.nDivisions_i && j< obj.nDivisions_j
                        queryIndeces = {[i,j],[i-1,j],[i,j+1],[i-1,j+1]};
                        for k = 1:4
                            currentVertexIndeces = queryIndeces{k};
                            verteces{k} = [obj.X_data(currentVertexIndeces(2),currentVertexIndeces(1)),...
                                           obj.Y_data(currentVertexIndeces(2),currentVertexIndeces(1)),...
                                           obj.Z_data(currentVertexIndeces(2),currentVertexIndeces(1))];
                        end

                        directionVector_i = verteces{2} - verteces{1};
                        directionVector_j = verteces{1} - verteces{3};
                        
                    elseif i < obj.nDivisions_i && j == obj.nDivisions_j
                        queryIndeces = {[i,j],[i+1,j],[i,j-1],[i+1,j-1]};
                        for k = 1:4
                            currentVertexIndeces = queryIndeces{k};
                            verteces{k} = [obj.X_data(currentVertexIndeces(2),currentVertexIndeces(1)),...
                                           obj.Y_data(currentVertexIndeces(2),currentVertexIndeces(1)),...
                                           obj.Z_data(currentVertexIndeces(2),currentVertexIndeces(1))];
                        end

                        directionVector_i = verteces{2} - verteces{1};
                        directionVector_j = verteces{3} - verteces{1};

                    else % i == obj.nDivisions_i && j == obj.nDivisions_j
                        queryIndeces = {[i,j],[i-1,j],[i,j-1],[i-1,j-1]};
                        for k = 1:4
                            currentVertexIndeces = queryIndeces{k};
                            verteces{k} = [obj.X_data(currentVertexIndeces(2),currentVertexIndeces(1)),...
                                           obj.Y_data(currentVertexIndeces(2),currentVertexIndeces(1)),...
                                           obj.Z_data(currentVertexIndeces(2),currentVertexIndeces(1))];
                        end

                        directionVector_i = verteces{2} - verteces{1};
                        directionVector_j = verteces{3} - verteces{1};
                    end

                    obj.dS_vec{j,i} = cross(directionVector_i,directionVector_j);
                    if any(isnan(obj.dS_vec{j,i}))
                        obj.dS_vec{j,i} = [0,0,0];
                    end
                end
            end            
        end

        function extrapolateNullLengthNormals(obj)
            for i = 1:obj.nDivisions_i
                for j=1:obj.nDivisions_j
                    if norm(obj.dS_vec{j,i}) == 0

                        interpDistance = 2;
                        dS_matrix = zeros(interpDistance,4);
                        if i == 1
                            for k = 1:interpDistance
                                dS_matrix(k,1) = k;
                                dS_matrix(k,[2,3,4]) = obj.dS_vec{j,i+k};
                            end
                        elseif i == obj.nDivisions_i
                            for k = 1:interpDistance
                                dS_matrix(k,1) = -k;
                                dS_matrix(k,[2,3,4]) = obj.dS_vec{j,i-k};
                            end
                        else
                            for k = 1:interpDistance/2
                                dS_matrix(k,1) = -k;
                                dS_matrix(k,[2,3,4]) = obj.dS_vec{j,i-k};
                            end
                            for k = interpDistance/2:interpDistance
                                dS_matrix(k,1) = k;
                                dS_matrix(k,[2,3,4]) = obj.dS_vec{j,i+k};
                            end 
                        end

                        dS_matrix_x = dS_matrix(:,2);
                        dS_matrix_y = dS_matrix(:,3);
                        dS_matrix_z = dS_matrix(:,4);

                        dS_vec_x = interp1(dS_matrix(:,1),dS_matrix_x,0,"linear","extrap");
                        dS_vec_y = interp1(dS_matrix(:,1),dS_matrix_y,0,"linear","extrap");
                        dS_vec_z = interp1(dS_matrix(:,1),dS_matrix_z,0,"linear","extrap");

                        obj.dS_vec{j,i} = [dS_vec_x,dS_vec_y,dS_vec_z];

                        if norm(obj.dS_vec{j,i}) == 0
                            %If norm is still 0
                            X_str = num2str(obj.X_data(j,i));
                            Y_str = num2str(obj.Y_data(j,i));
                            Z_str = num2str(obj.Z_data(j,i));
                            warning(['Singular point detected in surface: ',obj.name, ' in position: [',X_str,', ',Y_str,', ', Z_str,']', ...
                                '. Results at this point might not be correct'])                            
                            if i == obj.nDivisions_i
                                obj.dS_vec{j,i} = obj.dS_vec{j,i-1};
                            else
                                obj.dS_vec{j,i} = obj.dS_vec{j,i+1};
                            end

                            %If after all this, norm is still 0.
                            %Singularity spans more than a point and should
                            %not be considered a good surface
                            if norm(obj.dS_vec{j,i}) == 0
                                error(['More than one singular point detected in surface: ',obj.name, ' around point: [',X_str,', ',Y_str,', ', Z_str,']', ...
                                '. Surface definition should be reconsidered'])
                            end
                        end
                    end
                end
            end
        end

        function [vector3] = getFace_dS(obj,j,i)
            vector3 = [0,0,0];

            if i == obj.nDivisions_i || j == obj.nDivisions_j
                vector3 = obj.dS_vec{j,i};
                return
            end
            dS_matrix = [obj.dS_vec{j,i};obj.dS_vec{j+1,i};obj.dS_vec{j,i+1};obj.dS_vec{j+1,i+1}];
            dSx = dS_matrix(:,1);
            dSy = dS_matrix(:,2);
            dSz = dS_matrix(:,3);

            vector3 = [mean(dSx),mean(dSy),mean(dSz)];
        end

        function CPface = getFace_CP(obj,j,i)
            if i == obj.nDivisions_i || j == obj.nDivisions_j
                CPface = obj.Cp_data(j,i);
                return
            end
            CPface = mean([obj.Cp_data(j,i),obj.Cp_data(j+1,i),obj.Cp_data(j,i+1),obj.Cp_data(j+1,i+1)]);
        end
    end
end