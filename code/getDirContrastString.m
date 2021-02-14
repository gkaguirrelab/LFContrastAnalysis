function dirConStr = getDirContrastString(direction,contrast)
% this function takes in the direction and contrast numerical coding used
% in the experiment for bookkeeping purposes and returns a string with the
% corresponding direction and contrast. 

switch num2str(direction)
    case ('1');
        dirStr = '-45 deg';
    case ('2');
        dirStr = '45 deg';
    case ('3');
        dirStr = '0 deg';
    case ('4');
        dirStr = '90 deg';
    case ('5');
        dirStr = '-22.5 deg';
    case ('6');
        dirStr = '22.5 deg';
    case ('7');
        dirStr = '67.5 deg';
    case ('8');
        dirStr = '112.5 deg';
end

switch num2str(contrast)
    case ('1');
        dirConStr = [dirStr  ' 100%'];
    case ('2');
        dirConStr = [dirStr  ' 50%'];
    case ('3');
        dirConStr = [dirStr  ' 25%'];
    case ('4');
        dirConStr = [dirStr  ' 12.5%'];
    case ('5');
        dirConStr = [dirStr  ' 6.125%'];
end