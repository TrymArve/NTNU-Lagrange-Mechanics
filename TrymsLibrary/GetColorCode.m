function[code] = GetColorCode(color,brightness)

% color      - single-quote string containing the name of a color, or a one
%              letter code
% brightness - Multiplies the color code by this factor(caps at 1)

CODES = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560], ...
          [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0.5 0 0.5], ...
          [1 1 1]*0.5,[0.9686 0.498 0.7451],[0 0 0],[1 1 1]};
letter = ['b','o','y','p','g','c','r','m','e','i','k','w'];
word   = ["blue","orange","yellow","purple","green","cyan","red","magenta","grey","pink","black","white"];

try

if (ischar(color)) || (isstring(color))
    color = char(color);
    if length(color) == 1
        ind = find(letter == color,1);
    else
        ind = find(word == color,1);
    end
else
    letter = ['b','o','y','p','g','c','r','m','e','i','k']; %Do not include white if cycling through colors
    ind = mod(color-1,length(letter))+1;
end

code = CODES{ind};

if nargin > 1
    code = code*brightness;
    ind = find(code > 1);
    code(ind) = 1;
end

catch 
    disp('ERROR using GetColorCode.')
    disp('Available colors:')
    disp(word')
    disp('single letter codes:')
    disp(letter')
    
end