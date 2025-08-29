function linescanGUI(linescan, filtered_peaks, img)
    % Create the figure window
    fig = figure('Name', 'Linescan and Image Viewer', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    
    % Create an axes for the linescan plot
    axLinescan = axes('Parent', fig, 'Position', [0.05, 0.1, 0.4, 0.8]);
    imagesc(axLinescan, linescan);
    xlabel(axLinescan, 'm/z values');
    ylabel(axLinescan, 'Position along feature');
    title(axLinescan, 'Linescan');
    colormap(axLinescan, 'hot');
    colorbar(axLinescan);
    
    % Create an axes for the full image
    axFullImage = axes('Parent', fig, 'Position', [0.55, 0.1, 0.4, 0.8]);

    % Create a text box to display the m/z value
    txtMZ = uicontrol('Style', 'text', 'Parent', fig, 'Position', [400, 550, 200, 30], 'FontSize', 12, 'String', 'm/z Value: ');

    % Button to save the image
    btnSave = uicontrol('Style', 'pushbutton', 'Parent', fig, 'Position', [400, 50, 200, 30], 'String', 'Save Image', 'Callback', @saveImage);

    % Set up data cursor mode to enable column selection in the linescan
    dcm = datacursormode(fig);
    datacursormode on;
    set(dcm, 'UpdateFcn', @cursorCallback);
    
    % Global variable to store the currently selected column
    selectedColumn = 1;
    
    % Callback function when a point is selected in the linescan
    function txt = cursorCallback(~, event_obj)
        % Get the column index from the linescan plot
        pos = get(event_obj, 'Position');
        selectedColumn = pos(1);  % Column corresponds to the m/z index
        
        % Display the m/z value in the text box
        mzValue = filtered_peaks(selectedColumn);
        set(txtMZ, 'String', sprintf('m/z Value: %.4f', mzValue));
        
        % Display the corresponding full image in the axFullImage axes
        fullImage = squeeze(img(:, :, selectedColumn));
        imagesc(axFullImage, fullImage);
        axis(axFullImage, 'image');
        colormap(axFullImage, 'hot');
        colorbar(axFullImage);
        title(axFullImage, sprintf('m/z = %.4f', mzValue));
        
        % Return the string for the data cursor
        txt = {['Column: ', num2str(selectedColumn)], ['m/z: ', num2str(mzValue)]};
    end

    % Callback function for saving the image
    function saveImage(~, ~)
        % Get the current m/z value and full image
        mzValue = filtered_peaks(selectedColumn);
        fullImage = squeeze(img(:, :, selectedColumn));
        
        % Open the save file dialog
        [fileName, pathName] = uiputfile('*.png', 'Save Image As');
        if fileName == 0
            return; % If the user cancels, do nothing
        end
        
        % Save the image as a PNG file with the appropriate title and colorbar
        figSave = figure('Visible', 'off'); % Create a hidden figure to save
        imagesc(fullImage);
        axis image;
        colormap('hot');
        colorbar;
        title(sprintf('Full Image for m/z = %.4f', mzValue));
        
        % Save the figure
        saveas(figSave, fullfile(pathName, fileName));
        close(figSave); % Close the temporary figure after saving
    end
end
