function [X,Y] = f_meshgrid(xVector,fun,yDivisions,options)

arguments
    xVector (1,:) {mustBeNumeric}
    fun function_handle
    yDivisions (1,1) {mustBeNumeric}
    options.Simetry {mustBeMember(options.Simetry,{'true','false'})} = 'false'
    options.functionAxis {mustBeMember(options.functionAxis,{'X','Y'})} = 'X'
    %options.computationMethod {mustBeMember(options.computationMethod,{'serial','parallel','gpuArray'})} = 'serial'
end

lengthX = length(xVector);
X = zeros(yDivisions,lengthX);
Y = zeros(yDivisions,lengthX);
tic

if options.Simetry == 0
    scaleFactor = linspace(0,1,yDivisions);
else
    scaleFactor = linspace(-0.5,0.5,yDivisions);
end

for i=1:lengthX
X(:,i) = xVector(i);
    for j=1:yDivisions
        Y(j,i) = scaleFactor(j)*fun(xVector(i));
    end
end

if strcmp(options.functionAxis,'Y')
    %Swap the results since now Y is the main axis
    [Y,X] = deal(X,Y);
end
t = toc;
disp(['2D mesh generated in: ',num2str(t*1000),' miliseconds'])
end
