% ToolContactArea(X,Y,Z) 
%
% Writes comma delimited triplets of points lying on contact area to a
% file specified via GUI. GUI launches from current working directory.
%
% Companion subfunction to AFFC.m
% M.J. Roy, 2016
function ToolContactPointCloudWrite(X_ToolContact,...
    Y_ToolContact,Z_ToolContact)


q=0;
realvalue = isfinite(X_ToolContact);
[m,n]=size(realvalue);
for j=1:m
    for k=1:n
        if realvalue(j,k) == 1
            q=q+1;
            PointsOutput(q,1)=X_ToolContact(j,k);
            PointsOutput(q,2)=Y_ToolContact(j,k);
            PointsOutput(q,3)=Z_ToolContact(j,k);
        end
    end
end
[FileName,PathName] = uiputfile('*.asc','Save point cloud as:',pwd);
dlmwrite(fullfile(PathName, FileName),PointsOutput);

end %function