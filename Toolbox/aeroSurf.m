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
        dSx;
        dSy;
        dSz;
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

            disp(['Normals calculated for surface: ', obj.name])
        end

        function flipNormals(obj)
            obj.dSx = transpose(obj.dSx);
            obj.dSy = transpose(obj.dSy);
            obj.dSz = transpose(obj.dSz);
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
            quiver3(obj.X_data,obj.Y_data,obj.Z_data,obj.dSx,obj.dSy,obj.dSz);
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
            [obj.dSx,obj.dSy,obj.dSz] = surfnorm(obj.X_data,obj.Y_data,obj.Z_data);
        end

        function [vector3] = getFace_dS(obj,j,i)
            vector3 = [0,0,0];
            dS_matrix = zeros(3,4);

            if i == obj.nDivisions_i || j == obj.nDivisions_j
                vector3 = [obj.dSx(j,i),obj.dSy(j,i),obj.dSz(j,i)];
                return
            end
            dS_matrix(1,:) = [obj.dSx(j,i),obj.dSx(j,i+1),obj.dSx(j+1,i),obj.dSx(j+1,i+1)];
            dS_matrix(2,:) = [obj.dSy(j,i),obj.dSy(j,i+1),obj.dSy(j+1,i),obj.dSy(j+1,i+1)];
            dS_matrix(3,:) = [obj.dSz(j,i),obj.dSz(j,i+1),obj.dSz(j+1,i),obj.dSz(j+1,i+1)];

            vector3 = [mean(dS_matrix(1,:)),mean(dS_matrix(2,:)),mean(dS_matrix(3,:))];
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