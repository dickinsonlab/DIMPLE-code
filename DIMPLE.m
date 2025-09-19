%Set data directory
[peaklist] = batchcdfread('your file',number of rows);
%Set threshold factor for enrichment over background intensity 
minimum_fold_enrichment = 1.5;
%Set hot pixel detection parameters, decreasing the threshold will make hot
%pixel suppression more aggressive
hot_pixel = false;
hot_thresh = 7.7;
%Set mass uncertainty for peak matching across pixels
mass_uncertainty = 0.0000025;


img_size = [length(peaklist),size(peaklist{1},1)];
sum_img = zeros(img_size);
spectra = cell(img_size);
for i = 1:img_size(1)
    for j = 1:img_size(2)
        spectrum = peaklist{i}{j};
        spectra{i,j} = spectrum;
        sum_img(i,j) = sum(spectrum(:,2));
    end
end

%Reduce spectra to 5 significant digits
for i = 1:img_size(1)
    for j = 1:img_size(2)
        spectrum = spectra{i, j};        
        % Round m/z values to 4 significant digits
        rounded_mz = round(spectrum(:, 1), 5);
        [mz, ~, idx] = unique(rounded_mz);  % Get unique m/z and their indices
        ints = accumarray(idx, spectrum(:, 2));  % Sum intensities for each unique m/z
        spectra{i, j} = [mz, ints];
    end
end

%hot pixel detection
if hot_pixel
    bw = imregionalmax(sum_img,8);
    bw = imclearborder(bw);
    bg_img = imerode(sum_img,strel('square',3));
    ratio = sum_img./bg_img;
    %Any pixel more than 3 times as bright as the dimmest pixel in its
    %8-neighborhood is flagged as a hot pixel
    hot = ratio>hot_thresh;
    hot = imclearborder(hot);
    [r,c] = find(hot);
    for i = 1:length(r)
        block = spectra(r(i)-1:r(i)+1,c(i)-1:c(i)+1);
        block = reshape(block,[1,9]);
        center = block{5};
        center_int = sum(center(:,2));
        block = block([1:4,6:9]);
        block_int = 0;
        for j = 1:8
            block_int = block_int + sum(block{j}(:,2));
        end
        block_int = block_int/8;
        %Scale the intensities of the center hot pixel such that the total
        %signal is equal to the mean total signal of the surrounding 8
        %pixels
        center(:,2) = center(:,2)/(center_int/block_int);
        spectra{r(i),c(i)} = center;
    end
    
    %Recompute the sum_img
    sum_img = zeros(img_size);
    for i = 1:img_size(1)
        for j = 1:img_size(2)
            spectrum = spectra{i,j};
            sum_img(i,j) = sum(spectrum(:,2));
        end
    end

end


%Run this block only once to create an initial mask

figure
imagesc(log10(sum_img))
axis image
colorbar
roi = drawpolygon;
mask = roi.createMask;
figure
imshow(mask)
[r,c] = find(mask);

%pseudo bulk sample spectrum
mz_all = [];
for i = 1:length(r)
    spectrum = spectra{r(i), c(i)};
    mz_all = [mz_all; spectrum(:, 1)];
end

mz = unique(mz_all, 'stable');
sum_spectrum = zeros(length(mz), 1); 
[~, mz_indices] = ismember(mz, mz); 

for i = 1:length(r)
    spectrum = spectra{r(i), c(i)};
    [~, idx_in_mz] = ismember(spectrum(:, 1), mz);
    sum_spectrum(idx_in_mz) = sum_spectrum(idx_in_mz) + spectrum(:, 2);
    
    fprintf("row %u col %u len %u\n", r(i), c(i), length(spectrum));
end
[mz,i] = sort(mz);
sum_spectrum = sum_spectrum(i);

filtered = medfilt1(sum_spectrum,501);
filtered_spectrum = sum_spectrum-filtered;
filtered_spectrum(filtered_spectrum<0)=0;
peaks = mspeaks(mz, filtered_spectrum,'MULTIPLIER',100,'SHOWPLOT',0);
peaks = peaks(~isnan(peaks(:,1)),:);
%Suppress peaks within 2.5 ppm prioritizing stronger peaks
[~,index] = sort(peaks(:,2),'descend');
for i = 1:length(index)
    center = peaks(index(i),1);
    uncertainty = center*mass_uncertainty;
    if center>0
        conflicts = peaks(:,1)>=(center - uncertainty) & peaks(:,1)<=(center + uncertainty);
        conflicts(index(i))=0;
        peaks(conflicts,1) = -1;
    end
end
sample = peaks(peaks(:,1)>0,:);
sample(:,2) = sample(:,2)./length(r);

[r, c] = find(~mask);

mz_all = [];
for i = 1:length(r)
    spectrum = spectra{r(i), c(i)};
    mz_all = [mz_all; spectrum(:, 1)];
end

mz = unique(mz_all, 'stable');
sum_spectrum = zeros(length(mz), 1);

for i = 1:length(r)
    spectrum = spectra{r(i), c(i)};
    [~, idx_in_mz] = ismember(spectrum(:, 1), mz);
    sum_spectrum(idx_in_mz) = sum_spectrum(idx_in_mz) + spectrum(:, 2);

    fprintf("row %u col %u len %u\n", r(i), c(i), length(spectrum));
end

[mz,i] = sort(mz);
sum_spectrum = sum_spectrum(i);

filtered = medfilt1(sum_spectrum,501);
filtered_spectrum = sum_spectrum-filtered;
filtered_spectrum(filtered_spectrum<0)=0;
peaks = mspeaks(mz, filtered_spectrum,'MULTIPLIER',100,'SHOWPLOT',0);
peaks = peaks(~isnan(peaks(:,1)),:);
% Sort peaks by intensity (descending)
[~, sorted_indices] = sort(peaks(:, 2), 'descend');

% Precompute uncertainties for each peak
uncertainties = peaks(:, 1) * mass_uncertainty;

for i = 1:length(sorted_indices)
    idx = sorted_indices(i); 
    center = peaks(idx, 1);

    if center < 0
        continue;
    end

    lower_bound = center - uncertainties(idx);
    upper_bound = center + uncertainties(idx);

    conflicts = (peaks(:, 1) >= lower_bound) & (peaks(:, 1) <= upper_bound);
    conflicts(idx) = false;
    peaks(conflicts, 1) = -1;
end
background = peaks(peaks(:,1)>0,:);
background(:,2) = background(:,2)./length(r);


%suppress sample peaks for enrichment over background
[m,~] = size(sample);
unique_peaks = sample;


%Suppress peaks below enrichment factor threshold:
for i = 1:m
    ind = find(sample(i,1)==background(:,1));
    if ~isempty(ind)
        if sample(i,2)/background(ind,2)<minimum_fold_enrichment
            unique_peaks(i,1) = -1;
        end
    end
end
unique_peaks = unique_peaks(unique_peaks(:,1)>0,:);

filtered_peaks = unique_peaks(:,1);
[r,c] = find(mask);

harmonized_spectra = cell(size(spectra));
ppm_tolerances = filtered_peaks * mass_uncertainty;

for i = 1:length(r)
    spectrum = spectra{r(i), c(i)};
    
    bw = spectrum(:, 2) > 0;
    [L, N] = bwlabel(bw);

    harm_spectrum = zeros(N, 2);
    peak_centers = zeros(N, 1);
    for j = 1:N
        peak_centers(j) = mean(spectrum(L == j, 1));
    end

    indices = knnsearch(filtered_peaks, peak_centers, "K", 1);

    % Collect valid peaks and their max intensities
    harm_spectrum_idx = 1; 
    for j = 1:N
        % Get the current peak and corresponding filtered_peak
        peak = peak_centers(j);
        ind = indices(j);
        filtered_peak = filtered_peaks(ind);
        
        % Check if the peak is within the 2.5 ppm range of the filtered_peak
        lower_bound = filtered_peak - ppm_tolerances(ind);
        upper_bound = filtered_peak + ppm_tolerances(ind);
        
        if peak >= lower_bound && peak <= upper_bound
            % Peak is valid, include it in the harmonized spectrum
            harm_spectrum(harm_spectrum_idx, :) = [filtered_peak, max(spectrum(L == j, 2))];
            harm_spectrum_idx = harm_spectrum_idx + 1;
        end
    end

    harm_spectrum = harm_spectrum(1:harm_spectrum_idx-1, :);
    full_spectrum = zeros(length(filtered_peaks), 1);
    [~, ia, ib] = intersect(filtered_peaks, harm_spectrum(:, 1));
    full_spectrum(ia) = harm_spectrum(ib, 2);
    harmonized_spectra{r(i), c(i)} = full_spectrum;

    i / length(r)
end


%Build the multispectral image
img = zeros([img_size,length(filtered_peaks)]);
for i = 1:length(r)
    img(r(i),c(i),:) = harmonized_spectra{r(i),c(i)};
end

%saves the 3D image as a floating point tiff stack for visualization in
%image viewers like Fiji / ImageJ
% for z = 1:length(final_peaks)
%     t = Tiff(strcat(num2str(z),'.tif'), 'w');
%     tagstruct.ImageLength = size(filtered_img(:,:,z), 1);
%     tagstruct.ImageWidth = size(filtered_img(:,:,z), 2);
%     tagstruct.Compression = Tiff.Compression.None;
%     tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
%     tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%     tagstruct.BitsPerSample = 32;
%     tagstruct.SamplesPerPixel = 1;
%     tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%     t.setTag(tagstruct);
%     t.write(single(filtered_img(:,:,z)));
%     t.close();
% end

%Many peaks show no structured localization, appearing in a salt and pepper
%pattern, this step removes these peaks from the set
valid = zeros(1,length(filtered_peaks));
for i = 1:length(filtered_peaks)
    im = imerode(img(:,:,i),strel('disk',1));
    if sum(im(:))>0
        valid(i) = 1;
    end
end

final_peaks = filtered_peaks(valid>0);
filtered_img = img(:,:,valid>0);


%linescan assuming root is straight and vertical
%if root is horizontal, sum across dimension 1
%if root is vertical, sum across dimension 2
linescan = squeeze(max(filtered_img,[],2));
[m,n] = size(linescan);
for i = 1:n
    %Smooth and normalize the linescans
    linescan(:,i) = medfilt1(linescan(:,i)./max(linescan(:,i)),10);
end

D = pdist(linescan');
Z = linkage(squareform(D),'ward');
dendrogram(Z)
T = cluster(Z,'maxclust',4);
[B,I] = sort(T);
figure;imagesc(linescan(:,I))

%Determine cluster assignment of peak
idx = knnsearch(final_peaks,179.0105,'K',1);
figure;imagesc(filtered_img(:,:,idx))
highlight = T(I(idx));

%Plot mean linscans for each cluster normalized to area under curve
figure;
hold on
labels = string(zeros(1,max(T)));
for i = 1:max(T)
    if i == highlight
        plot(median(linescan(:,T==i)')./sum(median(linescan(:,T==i)'),'all'),'LineWidth',3);
    else
        plot(median(linescan(:,T==i)')./sum(median(linescan(:,T==i)'),'all'),'LineWidth',1.5);
    end
    labels(i) = strcat('Cluster ',num2str(i));
end
legend(labels);
hold off

%To have linescan annotated with cluster colormap
 figure;

% Create first axes for cluster assignments
ax1 = axes('Position', [0.1, 0.75, 0.8, 0.15]); % [left, bottom, width, height]
imagesc(B');
colormap(ax1, cool); % Set a distinct colormap just for this subplot
yticks([]);
set(ax1, 'XTick', []);
title('Cluster Assignments');

% Create second axes for the linescan image
ax2 = axes('Position', [0.1, 0.1, 0.8, 0.6]);
imagesc(linescan(:, I));
colormap(ax2, parula); % Set a different colormap here
xlabel('Sorted Indices');
ylabel('Line Position');
title('Sorted Linescan');

%To have dendrogram with cluster assignments aligned

figure;

% ===== 1. Dendrogram =====
ax3 = axes('Position', [0.1, 0.3, 0.8, 0.6]);
[H, ~, perm] = dendrogram(Z, 0);
set(ax3, 'XTick', []);
title('Dendrogram');

% ===== 2. Cluster Assignment Row =====
ax4 = axes('Position', [0.1, 0.2, 0.8, 0.05]);
imagesc(T(perm)');
colormap(ax4, jet); % or any distinct colormap
set(ax4, 'XTick', [], 'YTick', []);
title('Cluster Assignments');


%Determine cluster assignment of peak
idx = knnsearch(final_peaks,179.0105,'K',1);
figure;imagesc(filtered_img(:,:,idx))
highlight = T(I(idx));

%Plot mean linescans for each cluster normalized to area under curve
figure;
hold on
labels = string(zeros(1,max(T)));
for i = 1:max(T)
    if i == highlight
        plot(median(linescan(:,T==i)')./sum(median(linescan(:,T==i)'),'all'),'LineWidth',3);
    else
        plot(median(linescan(:,T==i)')./sum(median(linescan(:,T==i)'),'all'),'LineWidth',1.5);
    end
    labels(i) = strcat('Cluster ',num2str(i));
end
legend(labels);
hold off


%Example search for peak of interest

idx = knnsearch(final_peaks,117.0193,'K',1);
figure; imagesc(filtered_img(:,:,idx))
title(strcat('m/z = ',num2str(final_peaks(idx))))
axis image
axis off
colormap hot


%Run the linescan explorer GUI
linescanGUI(linescan(:,I),final_peaks(I),filtered_img(:,:,I));
