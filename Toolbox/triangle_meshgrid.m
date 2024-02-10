function [X,Y] = triangle_meshgrid(xVector,angleVector)
lengthX = length(xVector);
lengthAngle = length(angleVector);

X = zeros(lengthAngle,lengthX);
Y = zeros(lengthAngle,lengthX);

for i=1:lengthX
    X(:,i) = xVector(i);
    for j=1:lengthAngle
        Y(j,i) = xVector(i)*tan(angleVector(j));
    end
end

end