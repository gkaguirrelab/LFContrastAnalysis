function contrast = vectorAngle2LMScontrast(angles,plane)

for ii = 1:length(angles)
    switch plane
        case 'LM'
            if angles(ii) == -45
                contrast(1,ii) = 1;
                contrast(2,ii) = -1;
            elseif angles(ii) == 45
                contrast(1,ii) = 1;
                contrast(2,ii) = 1;
            else
                contrast(1,ii) = cosd(angles(ii));
                contrast(2,ii) = sind(angles(ii));
            end
            contrast(3,ii) = 0;
        case 'MS'
            % this needs to be though about more. only need LM plane now.
        case 'LS'
            % this needs to be though about more. only need LM plane now.
    end
end
end