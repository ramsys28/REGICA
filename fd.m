function [D] = fd(wave, method, param1,param2)
% Calculates fractal dimension of a waveform using either Katz's [1] or
% Sevcik's [2] algorithm.
%
% References:
% [1] Katz, M.J., Fractals and the analysis of waveforms, 
% Comput.Biol.Med. 18: 145, 1988
% [2] Sevcik, C., A procedure to Estimate the Fractal Dimension of
% Waveforms, Complexity International, volume 5, 1998, Available online:
% http://journal-ci.csse.monash.edu.au/ci/vol05/sevcik/

% Copyright (C) <2005>  <German Gomez-Herrero>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

TOL = 1e-6;

if nargin < 2, method = 'sevcik'; end


switch(lower(method)),
    case 'katz'
        n = length(wave);
        x = 1:n;
        y = wave;
        % Calculate the diameter
        d = sqrt((x-x(1)).^2+(y-y(1)).^2);
        d = max(d);
        % Calculate the length of the wave
        x = [ones(1,(n-1))];
        y = [wave(2:n)-wave(1:(n-1))];
        L = sum(sqrt(x.^2+y.^2));
        D = log10(n)/(log10(d/L)+log10(n));

    case 'sevcik',
        n = length(wave);
        x = 1:n;
        y = wave;
        % Map the wave to the unit square throught a double linear
        % transformation   
        span = (max(y)-min(y));
        if span < TOL,
            D = 1;
            return; 
        end
        y = (y-max(y))./span;         

        % calculate the length of the wave
        x = [(1/(n-1))*ones(1,(n-1))];
        y = [y(2:n)-y(1:(n-1))];
        L = sum(sqrt(x.^2+y.^2));
        D = 1+log(L)/log(2*(n-1));
        
    case 'sevcik_var',
        N = length(wave);
        wl = param1;
        ws = param2;
        ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        for i = 1:ne
           D(i) = fd(wave(init(i):min(final(i),N)),'sevcik'); 
        end
        D = var(D);
        
    case 'katz_var',
        N = length(wave);
        wl = param1;
        ws = param2;
        ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        for i = 1:ne
            D(i) = fd(wave(init(i):min(final(i),N)),'katz');
        end
        D = var(D);
        
        
    case 'sevcik_mean',
        N = length(wave);
        wl = param1;
        ws = param2;
        ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        for i = 1:ne
            D(i) = fd(wave(init(i):min(final(i),N)),'sevcik');
        end
        D = mean(D);  
        
    case 'katz_mean',
        N = length(wave);
        wl = param1;
        ws = param2;
        ovlength = wl-ws;
        init = 1:ws:length(wave);
        final = init+wl-1;
        ne = length(init);
        for i = 1:ne
            D(i) = fd(wave(init(i):min(final(i),N)),'katz');
        end
        D = mean(D);


    otherwise
        error('Unknown method');


end

