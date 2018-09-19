function contrast = vectorAngle2LMScontrast(angles,plane)

for ii = 1:length(angles)
    switch plane
        case 'LM'
            contrast(1,ii) = cosd(angles(ii));
            contrast(2,ii) = sind(angles(ii));
            contrast(3,ii) = 0;
        case 'MS'
            % this needs to be though about more. only need LM plane now.
        case 'LS'
            % this needs to be though about more. only need LM plane now.
    end
end
end