bodyLength = 5;

xdiv = 100;
ydiv = 100;

x = linspace(0,bodyLength,xdiv);
y = @(x) acosh(x+1)/2;

[X,Y] = f_meshgrid(x,y,ydiv,'Simetry','true');
Y = 2*Y;

lowerSeamHeight = @(x) lowerSeam(x,bodyLength);

Z_lower = zeros(size(X));
for i = 1:xdiv
    for j = 1:ydiv
        if Y(j,i) <= 0
            Z_lower(j,i) = (lowerSeamHeight(X(j,i))/sqrt(y(X(j,i))))*sqrt(Y(j,i)+y(X(j,i)));
        else
            Z_lower(j,i) = (lowerSeamHeight(X(j,i))/sqrt(y(X(j,i))))*sqrt(-Y(j,i)+y(X(j,i)));
        end
    end
end

upperSeamHeight = @(x) upperSeam(x,bodyLength);
interpFun = @(x) x^5;

Z_upper = zeros(size(X));
for i = 1:xdiv
    for j = 1:ydiv
        if Y(j,i) <= 0
            fun1 = ((upperSeamHeight(X(j,i)))/(y(X(j,i))))*Y(j,i)+upperSeamHeight(X(j,i));
            fun2 = ((Y(j,i)+y(X(j,i))).^4)/(y(X(j,i)).^4)*upperSeamHeight(X(j,i));
        else
            fun1 = ((-upperSeamHeight(X(j,i)))/(y(X(j,i))))*Y(j,i)+upperSeamHeight(X(j,i));
            fun2 = ((Y(j,i)-y(X(j,i))).^4)/(y(X(j,i)).^4)*upperSeamHeight(X(j,i));
        end
        Z_upper(j,i) = fun1 * (1-interpFun(X(j,i)/bodyLength)) + fun2 * (interpFun(X(j,i)/bodyLength)); 
    end
end

Z_lower(isnan(Z_lower)) = 0;
Z_upper(isnan(Z_upper)) = 0;
HipersonicBody = aeroShape(X,Y,Z_lower,'lower',X,Y,Z_upper,'upper');

HipersonicBody.flipNormals('lower')
%HipersonicBody.drawNormals

aoa = deg2rad(0);
beta = deg2rad(60);

tic
HipersonicNewton_solveCP(HipersonicBody,aoa,beta,'Cpt2',1.92)
HipersonicBody.drawCP
CA =[0.25*bodyLength,0,0];
[CF,CM] = HipersonicBody.solveForcesAndMoments(CA)
toc

%% Solve polars

n_points = 100;
angleRecorrido = linspace(-deg2rad(90),deg2rad(90),n_points);
CLrecorrido = zeros(n_points,1);
CDrecorrido = zeros(n_points,1);
CMalpharecorrido = zeros(n_points,1);


beta = 0;

tic
parfor i = 1:length(angleRecorrido)
    aoa = angleRecorrido(i);

    HipersonicNewton_solveCP(HipersonicBody,aoa,beta);
    [CF,CM] = HipersonicBody.solveForcesAndMoments(CA);

    Rotation_XYZtoVel = [cos(aoa)*cos(beta),cos(aoa)*sin(beta),sin(aoa);
                         -sin(beta),          cos(beta),           0;
                         -sin(aoa)*cos(beta), -sin(aoa)*sin(beta), cos(aoa)];

    CF_vel = Rotation_XYZtoVel * CF;
    CM_vel = Rotation_XYZtoVel * CM;
    CLrecorrido(i) = CF_vel(3);
    CDrecorrido(i) = CF_vel(1);
    CMalpharecorrido(i) = CM_vel(2);
end
toc

%% PLOTS
figure('Name','CL')
plot(rad2deg(angleRecorrido),CLrecorrido,'DisplayName','CL')
grid
legend

figure('Name','CD')
plot(rad2deg(angleRecorrido),CDrecorrido,'DisplayName','CD')
grid
legend

figure('Name','CL/CD')
plot(rad2deg(angleRecorrido),CLrecorrido./CDrecorrido,'DisplayName','CL/CD')
grid
legend

figure('Name','CMalpha')
plot(rad2deg(angleRecorrido),CMalpharecorrido,'DisplayName','CM_\alpha')
grid
legend
%% Piecewise functions
function z = lowerSeam(x,length)
    if x < (3/5)*length
        z = -0.55*sqrt(x);
    else
        z = (-0.55*sqrt((3/5)*length)/((3/5)*length))*x;
    end
end

function z = upperSeam(x,length)
    if x < (4/5)*length
        z = (0.5/length)*x;
    else
        z = (1.75/length)*x-1;
    end
end

