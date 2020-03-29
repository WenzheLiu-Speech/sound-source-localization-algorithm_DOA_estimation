function [pfEstAngles,figHandle] = post_findPeaks(specGlobal, azimuth, elevation, azimuthGrid, elevationGrid, nsrc, minAngle, displayResults)
if(displayResults)
    figHandle = figure;
else
    figHandle = -1;
end
if(isempty(minAngle))
    minAngle = 1;
elseif(minAngle < abs(azimuthGrid(2)- azimuthGrid(1)) && nsrc>1)
    error('Error[findPeaks]:谱峰之间最小夹角必须大于方位角/俯仰角的分辨率');
else
    
end

nAzi = length(azimuth);
nEl = length(elevation);

if(nAzi == 1 || nEl == 1) % 水平面或竖直面
    [~, ind] = findpeaks(specGlobal, 'minpeakdistance',minAngle, 'sortstr','descend');
    numberSourcesFound = min(nsrc,length(ind));
    
    if(nAzi == 1)
        pfEstAngles = [azimuth.*ones(numberSourcesFound,1) elevationGrid(ind(1:numberSourcesFound))'];
    else
        pfEstAngles = [azimuthGrid(ind(1:numberSourcesFound))' elevation.*ones(numberSourcesFound,1)];
    end
    
    % Display
    if (displayResults)
        figure(figHandle);
        if(nAzi == 1)
            plot(elevationGrid,specGlobal);
            hold on
            plot(elevationGrid(ind(1:numberSourcesFound)),specGlobal(ind(1:numberSourcesFound)),'*k','MarkerSize',15,'linewidth',1.5);
            xlabel('\phi (degrees)');
            ylabel('Angular spectrum');
            hold off;
            title(['\Phi(\theta,\phi) for \theta = ' num2str(azimuth) '  |  markers : sources found ']);
        else
            plot(azimuthGrid,specGlobal);
            hold on
            %display sources found
            plot(azimuthGrid(ind(1:numberSourcesFound)),specGlobal(ind(1:numberSourcesFound)),'*k','MarkerSize',15,'linewidth',1.5);
            xlabel('\theta (degrees)');
            ylabel('Angular spectrum');
            hold off;
            title(['\Phi(\theta,\phi) for \phi = ' num2str(elevation) '  |  markers : sources found ']);
        end
    end
    % 
else  % 二维
    % Convert angular spectrum in 2D
    ppfSpec2D = (reshape(specGlobal,nAzi,nEl))';
    ppfPadpeakFilter = ones(size(ppfSpec2D,1)+2,size(ppfSpec2D,2)+2) * -Inf;
    ppfPadpeakFilter(2:end-1,2:end-1) = ppfSpec2D;

    % Find peaks : compare values with their neighbours
    ppiPeaks = ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,2:end-1) & ... % top
        ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  2:end-1) & ... % bottom
        ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(2:end-1,1:end-2) & ... % right
        ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(2:end-1,3:end)   & ... % left
        ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,1:end-2) & ... % top/left
        ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,3:end)   & ... % top/right
        ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  1:end-2) & ... % bottom/left
        ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  3:end);        % bottom/right

    % number of local maxima
    iNbLocalmaxima = sum(sum(ppiPeaks));

    % local maxima with corrresponding values
    ppfSpec2D_peaks = (ppfSpec2D - min(min(ppfSpec2D))) .* ppiPeaks; % substract min value : avoid issues (when sorting peaks) if some peaks values are negatives

    % sort values of local maxima
    pfSpec1D_peaks= reshape(ppfSpec2D_peaks',1,nEl*nAzi);
    [~,piIndexPeaks1D] = sort(pfSpec1D_peaks,'descend');

    piEstSourcesIndex = piIndexPeaks1D(1);  % first source is the global maximum (first one in piSortedPeaksIndex1D)
    index = 2; % search index in piSortedPeaksIndex1D
    numberSourcesFound = 1; % set to one as global maximum is already selected as source

    %Filter the list of peaks found with respect to minAngle parameter
    while (numberSourcesFound < nsrc && index <= iNbLocalmaxima)

        bAngleAllowed = 1;
        % verify that current direction is allowed with respect to minAngle and sources already selected
        for i = 1:length(piEstSourcesIndex)

            % distance calculated using curvilinear abscissa (degrees) - ref. : http://geodesie.ign.fr/contenu/fichiers/Distance_longitude_latitude.pdf
            dist=acosd(sind(elevationGrid(piEstSourcesIndex(i)))*sind(elevationGrid(piIndexPeaks1D(index)))+cosd(elevationGrid(piEstSourcesIndex(i)))*cosd(elevationGrid(piIndexPeaks1D(index)))*cosd(azimuthGrid(piIndexPeaks1D(index))-azimuthGrid(piEstSourcesIndex(i))) );

            if(dist <minAngle)
                bAngleAllowed =0;
                break;
            end
        end

        % store new source
        if(bAngleAllowed)
            piEstSourcesIndex = [piEstSourcesIndex,piIndexPeaks1D(index)];
            numberSourcesFound = numberSourcesFound +1;
        end

        index = index + 1;
    end

    pfEstAngles = [azimuthGrid(piEstSourcesIndex)' elevationGrid(piEstSourcesIndex)'];

%% Display results
    if (displayResults)
        figure(figHandle);
        colormap(jet); % bleu jaune rouge
        imagesc(azimuth,elevation,ppfSpec2D);
        set(gca,'YDir','normal');
        hold on;

        %display sources found
        for i =1:length(pfEstAngles(:,1))
            handle=plot(pfEstAngles(i,1),pfEstAngles(i,2),'*k','MarkerSize',15,'linewidth',1.5);
        end

        xlabel('\theta (degrees)');
        ylabel('\phi (degrees)');
        hold off;
        title('\Phi(\theta,\phi)   |  markers : sources found ');
    end
end
drawnow;
end