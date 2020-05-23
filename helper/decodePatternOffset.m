function [ patternOffsetPosition ] = decodePatternOffset( patternNumber )
%DECODEPATTERNOFFSET solves for if there was an offset between x-position
%diremtion of the pattern and where the object, likely bar was actually
%located
%
%   INPUT
% patternNumber the value stored in stimulus.panelParams.patternNum;
%   OUTPUT
% patternOffsetPosition - value in LED positions that the bar was offset
%       from the pattern dimentions
%       ie if patternOffsetPosition = 20, then when the pattern is in
%       x-position 1, the bar was at x position 20 on the panels
%   3/29/18 Yvette Fisher

if( patternNumber == 3 || patternNumber == 37 || patternNumber == 38 || patternNumber == 14 || patternNumber == 26 ||  patternNumber == 29 || patternNumber == 32 || patternNumber == 45 )
    
    patternOffsetPosition = 0; %NO offset
    
elseif( patternNumber == 15 || patternNumber == 27 ||  patternNumber == 30 || patternNumber == 33 || patternNumber == 46)
    
    patternOffsetPosition = 24; % 90 deg offset
    
elseif( patternNumber == 16 || patternNumber == 28  ||  patternNumber == 31 || patternNumber == 34 || patternNumber == 47)
    
    patternOffsetPosition = 48; % 180 deg offset
    
elseif( patternNumber == 23 )
    
    patternOffsetPosition = 0; %NO offset
    %TODO change to account for the inverted gain
else
    patternOffsetPosition = 0; %NO offset
    warning('ERROR pattern number was not set up to be decoded, offset defaulted to = 0');
end

end

