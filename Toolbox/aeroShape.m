classdef aeroShape < handle

    properties
        surfacesArray (1,:) aeroSurf
        surfacesNames (1,:) cell
        numberOfSurfaces uint16
    end

    methods
        function obj = aeroShape(X,Y,Z,name)
            arguments (Repeating)
                X {mustBeNumeric}
                Y {mustBeNumeric}
                Z {mustBeNumeric}
                name {mustBeTextScalar}
            end

            %We check if names are unique. If not, the function throws an
            %error
            obj.assertUniqueNames(name)

            obj.numberOfSurfaces = length(X);
            for i = 1:obj.numberOfSurfaces
                obj.surfacesArray(i) = aeroSurf();
                obj.surfacesArray(i).setGeometry(X{i},Y{i},Z{i},name{i});
                obj.surfacesNames{i} = name{i};
            end
        end

        function setSurfaces(obj,X,Y,Z,name)
           arguments
               obj aeroShape
           end
           arguments (Repeating)
               X {mustBeNumeric}
               Y {mustBeNumeric}
               Z {mustBeNumeric}
               name {mustBeTextScalar}
           end

           %Reset the object
           obj.emptyObject

           %We check if names are unique. If not, the function throws an
           %error
           obj.assertUniqueNames(name)

           obj.numberOfSurfaces = length(X);
           for i = 1:obj.numberOfSurfaces
               obj.surfacesArray(i) = aeroSurf();
               obj.surfacesArray(i).setGeometry(X{i},Y{i},Z{i},name{i});
               obj.surfacesNames{i} = name{i};
           end
        end

        function addSurfaces(obj,X,Y,Z,name)
           arguments
               obj aeroShape
           end
           arguments (Repeating)
               X {mustBeNumeric}
               Y {mustBeNumeric}
               Z {mustBeNumeric}
               name {mustBeTextScalar}
           end

           if nargin < 2
               error('addSurfaces needs to add at least one surface')
           end

           %We check if names are unique. If not, the function throws an
           %error
           obj.assertUniqueNames(name)
           
           previousNumberOfSurfaces = obj.numberOfSurfaces;
           numberOfNewSurfaces = length(X);
           obj.numberOfSurfaces = previousNumberOfSurfaces + numberOfNewSurfaces;

           %Increase the size of the array
           for i = 1:numberOfNewSurfaces
               obj.surfacesArray(end+i) = aeroSurf();
           end

           %Populate array
           for i = 1:numberOfNewSurfaces
               obj.surfacesArray(previousNumberOfSurfaces+i).setGeometry(X{i},Y{i},Z{i},name{i});
               obj.surfacesNames{previousNumberOfSurfaces+i} = name{i};
           end
        end

        function [surface,index] = getSurface(obj,name)
            arguments
                obj aeroShape
                name {mustBeTextScalar}
            end

            index = find(strcmp(obj.surfacesNames,name));

            if isempty(index)
                error(['Surface: ', name, ' not found'])
            end

            surface = obj.surfacesArray(index);
        end

        function deleteSurface(obj,name)
           arguments
                obj aeroShape
                name {mustBeTextScalar}
           end

            [~,index] = obj.getSurface(name);
            obj.surfacesArray(index) = [];
            obj.surfacesNames(index) = [];
            obj.numberOfSurfaces = obj.numberOfSurfaces - 1;
        end

        function flipNormals(obj,name)
           arguments
                obj aeroShape
                name {mustBeTextScalar}
           end
           [~,index] = obj.getSurface(name);
           obj.surfacesArray(index).flipNormals;
        end

        function setPremadeShape(obj,options)
            arguments
                obj aeroShape
                options.premade {mustBeMember(options.premade,['Paraboloid']),mustBeNonempty}
                options.length {mustBeScalarOrEmpty,mustBeNumeric} = 1;
                options.nDivisions_X {mustBeScalarOrEmpty,mustBeNumeric} = 10;
                options.nDivisions_Y {mustBeScalarOrEmpty,mustBeNumeric} = 10;
            end

            if strcmp(options.premade,'Paraboloid')
                xVec = linspace(0,options.length,options.nDivisions_X);
                y = @(x) real(sqrt(x));
                [X,Y] = f_meshgrid(xVec,y,options.nDivisions_Y,'Simetry','true');
    
                r = y(X)/2;
                Z_upper = sqrt(r.^2-Y.^2);
                Z_lower = -sqrt(r.^2-Y.^2);

                obj.setSurfaces(X,Y,Z_upper,'upper',X,Y,Z_lower,'lower');
                obj.surfacesArray(2).flipNormals;
            end
        end

        function [CF] = solveForces(obj)
            CF = [0;0;0];
            for i = 1:obj.numberOfSurfaces
                if ~all(size(obj.surfacesArray(i).Cp_data))
                    error(['CP_data not populated for surface: ',obj.surfacesArray(i).name,'. Use obj.solveCP() to populate CP field'])
                end
                CF_surf = obj.surfacesArray(i).getForces;
                CF = CF + CF_surf;
            end
        end

        function [CF,CM] = solveForcesAndMoments(obj,momentsPoint)
            CF = [0;0;0];
            CM = [0;0;0];

            for i = 1:obj.numberOfSurfaces
                if ~all(size(obj.surfacesArray(i).Cp_data))
                    error(['CP_data not populated for surface: ',obj.surfacesArray(i).name,'. Use obj.solveCP() to populate CP field'])
                end
                [CF_surf,CM_surf] = obj.surfacesArray(i).getForcesAndMoments(momentsPoint);
                CF = CF + CF_surf;
                CM = CM + CM_surf;
            end
        end

        function drawSurfaces(obj)
            figure
            for i = 1:obj.numberOfSurfaces
                obj.surfacesArray(i).drawGeometry
            end
        end

        function drawNormals(obj)
            figure
            for i = 1:obj.numberOfSurfaces
                obj.surfacesArray(i).drawNormals
            end
        end

        function drawCP(obj)
            figure
            for i = 1:obj.numberOfSurfaces
                obj.surfacesArray(i).drawCP
            end
        end
    end

    methods (Access = protected)
        function assertUniqueNames(obj,name)
            arguments
                obj aeroShape
                name cell {mustBeText}
            end

            if length(unique(name))~=length(name)
                error('Surfaces names passed in are invalid. Names must be unique in each shape object')
            end
            if ~isempty(obj.surfacesNames)
                if any(strcmp(name,obj.surfacesNames))
                    error(['Surface name: ', name{1}, ' is invalid. Names must be unique in each shape object'])
                end
            end
        end


        function emptyObject(obj)
            obj.surfacesArray = aeroSurf.empty;
            obj.surfacesNames = {};
            obj.numberOfSurfaces = 0;
        end
    end
end