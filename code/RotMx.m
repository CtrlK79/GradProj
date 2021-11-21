% Generate rotational matrix
function m = RotMx(axis, degree)
    if axis == 1
        m = [1 0 0; 0 cosd(degree) -sind(degree); 0 sind(degree) cosd(degree)];
    elseif axis == 2
        m = [cosd(degree) 0 sind(degree); 0 1 0; -sind(degree) 0 cosd(degree)];
    elseif axis == 3
        m = [cosd(degree) -sind(degree) 0; sind(degree) cosd(degree) 0; 0 0 1];
    elseif axis == 4
        m = [cosd(degree) -sind(degree); sind(degree) cosd(degree)];
    end
end
            