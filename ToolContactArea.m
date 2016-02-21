% [Axyz Axy Axz Ayz] = ToolContactArea(X,Y,Z) 
%
% Returns areas calculated from X, Y and Z component matrices by
% triangulation.
%
% Companion subfunction to AFFC.m
% M.J. Roy, 2016
function [AreaXYZ, AreaXY, AreaXZ, AreaYZ]= ToolContactArea(...
    X_ToolContact,Y_ToolContact,Z_ToolContact)

triInd = makeind(X_ToolContact);

X = X_ToolContact(triInd);
Y = Y_ToolContact(triInd);
Z = Z_ToolContact(triInd);

AreaXYZ = 0;
for i = 1:length(X(:,1))
    if sum(isnan(X(i,:)))==0
        A = [[Y(i,1) Z(i,1) 1]; [Y(i,2) Z(i,2) 1]; [Y(i,3) Z(i,3) 1];];
        B = [[Z(i,1) X(i,1) 1]; [Z(i,2) X(i,2) 1]; [Z(i,3) X(i,3) 1];];
        C = [[X(i,1) Y(i,1) 1]; [X(i,2) Y(i,2) 1]; [X(i,3) Y(i,3) 1];];
        
        A_e = .5*sqrt(det(A)^2+det(B)^2+det(C)^2);
        
        AreaXYZ = AreaXYZ+A_e;
    end
end


X = X_ToolContact(triInd);
Y = Y_ToolContact(triInd);
Z = zeros(size(Z_ToolContact(triInd)));

AreaXY = 0;
for i = 1:length(X(:,1))
    if sum(isnan(X(i,:)))==0
        A = [[Y(i,1) Z(i,1) 1]; [Y(i,2) Z(i,2) 1]; [Y(i,3) Z(i,3) 1];];
        B = [[Z(i,1) X(i,1) 1]; [Z(i,2) X(i,2) 1]; [Z(i,3) X(i,3) 1];];
        C = [[X(i,1) Y(i,1) 1]; [X(i,2) Y(i,2) 1]; [X(i,3) Y(i,3) 1];];
        
        A_e = .5*sqrt(det(A)^2+det(B)^2+det(C)^2);
        
        AreaXY = AreaXY+A_e;
    end
end


X = X_ToolContact(triInd);
Y = zeros(size(Y_ToolContact(triInd)));
Z = Z_ToolContact(triInd);

AreaXZ = 0;
for i = 1:length(X(:,1))
    if sum(isnan(X(i,:)))==0
        A = [[Y(i,1) Z(i,1) 1]; [Y(i,2) Z(i,2) 1]; [Y(i,3) Z(i,3) 1];];
        B = [[Z(i,1) X(i,1) 1]; [Z(i,2) X(i,2) 1]; [Z(i,3) X(i,3) 1];];
        C = [[X(i,1) Y(i,1) 1]; [X(i,2) Y(i,2) 1]; [X(i,3) Y(i,3) 1];];
        
        A_e = .5*sqrt(det(A)^2+det(B)^2+det(C)^2);
        
        AreaXZ = AreaXZ+A_e;
    end
end


X = zeros(size(X_ToolContact(triInd)));
Y = Y_ToolContact(triInd);
Z = Z_ToolContact(triInd);

AreaYZ = 0;
for i = 1:length(Y(:,1))
    if sum(isnan(Y(i,:)))==0
        A = [[Y(i,1) Z(i,1) 1]; [Y(i,2) Z(i,2) 1]; [Y(i,3) Z(i,3) 1];];
        B = [[Z(i,1) X(i,1) 1]; [Z(i,2) X(i,2) 1]; [Z(i,3) X(i,3) 1];];
        C = [[X(i,1) Y(i,1) 1]; [X(i,2) Y(i,2) 1]; [X(i,3) Y(i,3) 1];];
        
        A_e = .5*sqrt(det(A)^2+det(B)^2+det(C)^2);
        
        AreaYZ = AreaYZ+A_e;
    end
end

%return areas to console
fprintf('\n');
fprintf('Calculated areas . . . \n')
fprintf('Overall\t\t|X plane\t|Y plane\t|Z plane\n');
fprintf('%f\t|%f\t|%f\t|%f\n',AreaXYZ,AreaYZ,AreaXZ,AreaXY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triInd=makeind(x)
[u,v] =size(x);
Tri=[];
for j=1:u-1
    tri_2 = [];
    for I=1:v-1
        Tri = [Tri; [sub2ind(size(x),j,I) sub2ind(size(x),j+1,I) ...
            sub2ind(size(x),j+1,I+1)]];
        tri_2 = [tri_2; [sub2ind(size(x),j,I) sub2ind(size(x),j+1,I+1) ...
            sub2ind(size(x),j,I+1)]];
    end
    Tri=[Tri; tri_2];
end
triInd=Tri;
end %makemesh



end %main function