function [windDecreaseFactorTotal,gridCoordinate_new,dataPlot] = wakeEffectDemo(OWFapprovedArea,ULML,VLML,windTurbine,resolution)
% this function is to generate a wake effect demo for drawing graphs
x0 = 0:resolution:OWFapprovedArea.length;
y0 = 0:resolution:OWFapprovedArea.height;
xWT = 0:OWFapprovedArea.lengthSquare:OWFapprovedArea.length;
yWT = 0:OWFapprovedArea.lengthSquare:OWFapprovedArea.height;
% shift
    % ULML: surface eastward wind (originates in the east and blows in a westward direction.); VLML: surface northward wind
theta = atan(VLML / ULML); % ↙
for iRow = 1:OWFapprovedArea.nRow
    for iCol = 1:OWFapprovedArea.nColumn
        WTcoordinate_new{iRow,iCol} = [cos(theta),-sin(theta);sin(theta),cos(theta)] * [xWT(iCol); yWT(iRow)];
    end
end
% calculate wake effect
windDecreaseFactorTotal = zeros(size(y0,2),size(x0,2));
gridCoordinate_new = cell(size(y0,2),size(x0,2));
for ix0 = 1:size(x0,2)
    for iy0 = 1:size(y0,2)
        % new grid coordinates
        gridCoordinate_new{iy0,ix0} = [cos(theta),-sin(theta);sin(theta),cos(theta)] * [x0(ix0); y0(iy0)];
        % wake effect from each TW
        windDecreaseFactor = zeros(OWFapprovedArea.nRow,OWFapprovedArea.nColumn);
        for iRow = 1:OWFapprovedArea.nRow
            for iCol = 1:OWFapprovedArea.nColumn
                delta_x = gridCoordinate_new{iy0,ix0}(1) - WTcoordinate_new{iRow,iCol}(1);
                delta_y = abs(gridCoordinate_new{iy0,ix0}(2) - WTcoordinate_new{iRow,iCol}(2));
                if delta_x < 0
                    windDecreaseFactor(iRow,iCol) = 0;
                else
                    sigma = (windTurbine.rotorRadius + 0.56 / (log(windTurbine.hubHeight) - log(windTurbine.roughness)) .* delta_x) / 2;
                    windDecreaseFactor(iRow,iCol) = (1 - sqrt(1 - windTurbine.C_T / 2 ./ (sigma / windTurbine.rotorRadius)^2) ) * exp(-delta_y.^2 / 2 / sigma.^2);
                end
            end
        end
        windDecreaseFactorTotal(iy0,ix0) = sum(sum(windDecreaseFactor.^2));
    end
end
counter = 0;
dataPlot = zeros(size(x0,2)*size(y0,2),3);
for ix0 = 1:size(x0,2)
    for iy0 = 1:size(y0,2)
        counter = counter + 1;
        dataPlot(counter,1) = gridCoordinate_new{iy0,ix0}(1); % 
        dataPlot(counter,2) = gridCoordinate_new{iy0,ix0}(2);
        dataPlot(counter,3) = windDecreaseFactorTotal(iy0,ix0);
    end
end