% =========================================================
% Write out coordinates for a Neuropixels metadata file.
% If there is an ~snsGeomMap, that will be used. If not, and the probe is
% supported, the geometry will be derived from the ~snsShankMap.
% outType = 3 will make a copy of the metadata file with the ~snsGeomMap
% appended.
%
% If you get a 'Probe not supported error', open your data in the
% SpikeGLX viewer (version 20230202 or later); and export a very short
% segment; the ~snsGeomMap will be appended to the metadata file created for
% the exported file. You can then use this script to make a KS channel map.
% 
% Format selected with the outType variable
%     0 for tab delimited text coordinate file: 
%       chan index, x in um, y in um, shank index 
%     1 for Kilosort or Kilosort2 channel map file;
%     2 for strings to paste into JRClust .prm file;
%     3 adds muxTbl, snsGeomMap and other metadata items to existing file
%       to be compatible with SpikeGLX 032623-phase 30 and newer
%     4 returns the coordinates as a matlab struct and does not write files
%
% File is saved in current working directory 
%
% Jennifer Colonell, Janelia Research Campus
 
function out=SGLXMetaToCoords(varargin)

% Output selection:
outType = 0;
bProcLF = 1; % only used with outType 3; append new fields to matching lf meta
out=[];

if (length(varargin) == 0)
    % Ask user for metadata file
    [metaName,path] = uigetfile('*.meta', 'Select AP Metadata File');
    % Plot shank view
    bPlot = 1;
else
    % calling from another script
    inputCell = varargin(1);
    metaFullPath = inputCell{1};   % full path to the input res file
    [path, tempName, ext] = fileparts(metaFullPath);
    metaName = sprintf('%s%s', tempName, ext);
    inputCell = varargin(2);
    bProcLF = inputCell{1};  % append new fields to matching lf meta?
    bPlot = 0;
end

% Parse input file to get the metadata structure
meta = ReadMeta(metaName, path);

% Get coordinates for saved channels from snsGeomMap, if present,
% otherwise from snsShankMap
if isfield(meta,'snsGeomMap')
    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta);
elseif isfield(meta,'snsShankMap')    
    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta);
end

if bPlot
    plotSaved(xCoord, yCoord, shankInd, meta);
end

% Build output name and write out file
[~,fname,~] = fileparts(metaName);

switch outType
    
    case 0      %tab delimited, chan, x, y, shank
        newName = [fname,'-siteCoords.txt'];
        fid = fopen( newName, 'w');
        nchan = numel(xCoord);
        for i = 1:nchan
            currX = shankInd(i)*shankPitch + xCoord(i);
            fprintf( fid, '%d\t%d\t%d\t%d\n', i-1, currX, yCoord(i), shankInd(i));
        end
        fclose(fid);
        
    case 1     %KS2 *.mat
        newName = [fname,'_kilosortChanMap.mat'];
        nchan = numel(xCoord);
        chanMap = (1:nchan)';
        chanMap0ind = chanMap - 1;
        connected = logical(connected);
        xcoords = shankInd*shankPitch + xCoord;   %KS2 not yet using kcoords, so x coord includes shank sep
        ycoords = yCoord; % variable names need to  match KS standard
        kcoords = shankInd + 1;     %KS1 uses kcoords to force templates to be on one shank
        name = fname;
        save( newName, 'chanMap', 'chanMap0ind', 'connected', 'name', 'xcoords', 'ycoords', 'kcoords' );
    
    case 2  %strings to copy into JRC prm file
       % for JRC, shankMap = array of shank indicies, 1 based 
       % siteLoc = array of (x,y) pairs
       % siteMap = order in file (channel index, 1-based) of the locations in siteLoc
       nchan = numel(xCoord); 
       chans = (1:nchan); 
       newName = [fname,'_forJRCprm.txt'];
       fid = fopen( newName, 'w' );

       fprintf( fid, 'shankMap = [' );       
       for i = 1:nchan-1
           fprintf( fid, '%d,', shankInd(i) + 1 ); 
       end
       fprintf( fid, '%d];\n',shankInd(nchan) + 1 );
       
       xCoord = shankInd*shankPitch + xCoord; 
       fprintf( fid, 'siteLoc = [' );
       for i = 1:nchan-1
           fprintf(fid, '%d,%d;', xCoord(i), yCoord(i));
       end
       fprintf( fid, '%d,%d];\n', xCoord(nchan), yCoord(nchan) );
       
       fprintf( fid, 'siteMap = [' );
       for i = 1:nchan-1
           fprintf( fid, '%d,', chans(i) );
       end
       fprintf( fid, '%d];\n', chans(nchan) );
       fclose(fid);

    case 3 % add metadata items, including muxTbl and snsGeom
        newName = [fname,'_orig.meta'];
        % rename the metadata file (preserves the modify date)
        movefile(fullfile(path,metaName), fullfile(path,newName));
        copyfile(fullfile(path,newName), fullfile(path,metaName));
        [imChan0apGainStr, imChan0lfGainStr, imAnyChanFullBandStr] = imroMetaItems(meta);
        muxTableStr = getMuxTable(meta);
        snsGeomStr = snsGeom(meta, shankInd, xCoord, yCoord, connected, nShank, shankPitch, shankWidth );
        % add items to the metadata file
        fid = fopen(fullfile(path,metaName), 'a');
        fprintf(fid,'%s\n',imChan0apGainStr);
        fprintf(fid,'%s\n',imChan0lfGainStr);
        fprintf(fid,'%s\n',imAnyChanFullBandStr);
        fprintf(fid,'%s\n',muxTableStr);
        fprintf(fid,'%s\n',snsGeomStr);
        fclose(fid);
        
        if bProcLF
            lf_fname = sprintf('%s.lf', extractBefore(fname,'.ap'));
            lf_metaName = sprintf('%s.meta',lf_fname);
            if isfile(fullfile(path,lf_metaName))
                newName = [lf_fname,'_orig.meta'];
                movefile(fullfile(path,lf_metaName), fullfile(path,newName));
                copyfile(fullfile(path,newName), fullfile(path,lf_metaName));
                % add items to the metadata file
                fid = fopen(fullfile(path,lf_metaName), 'a');
                fprintf(fid,'%s\n',imChan0apGainStr);
                fprintf(fid,'%s\n',imChan0lfGainStr);
                fprintf(fid,'%s\n',imAnyChanFullBandStr);
                fprintf(fid,'%s\n',muxTableStr);
                fprintf(fid,'%s\n',snsGeomStr);
                fclose(fid);
            end
        end
        
    case 4      %matlab struct
        out = struct('nShank',nShank,'shankWidth',shankWidth,'shankPitch',shankPitch,...
                     'shankInd',shankInd,'xCoord',xCoord,'yCoord',yCoord,'connected',connected);
                 
end
end % SGLXMetaToCoords


% =========================================================
% Read metadata file and return meta structure,
% each tag becomes a field.
%
function [meta] = ReadMeta(metaName, path)
    % Parse ini file into cell entries C{1}{i} = C{2}{i}
    fid = fopen(fullfile(path, metaName), 'r');
% -------------------------------------------------------------
%    Need 'BufSize' adjustment for MATLAB earlier than 2014
%    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
    C = textscan(fid, '%[^=] = %[^\r\n]');
% -------------------------------------------------------------
    fclose(fid);

    % New empty struct
    meta = struct();

    % Convert each cell entry into a struct entry
    for i = 1:length(C{1})
        tag = C{1}{i};
        if tag(1) == '~'
            % remake tag excluding first character
            tag = sprintf('%s', tag(2:end));
        end
        meta = setfield(meta, tag, C{2}{i});
    end
end % ReadMeta


% =========================================================
% Return counts of each imec channel type that compose
% the timepoints stored in binary file.
%
function [AP,LF,SY] = ChannelCountsIM(meta)
    M = str2num(meta.snsApLfSy);
    AP = M(1);
    LF = M(2);
    SY = M(3);
end % ChannelCountsIM


% =========================================================
% Plot x z positions of all electrodes and saved channels
%
function plotSaved(xCoord, yCoord, shankInd, meta)
    % get geometry
    g = getGeomParams(meta);

    % calculate positions on one shank
    nCol = g.elecPerShank/g.rowsPerShank;
    rowInd = 0:g.elecPerShank-1;
    rowInd = floor(rowInd./nCol);
    oddRows = logical(mod(rowInd,2));
    evenRows = ~oddRows;

    colInd = 0:g.elecPerShank-1;
    colInd = mod(colInd,nCol);

    xall = colInd*g.horzPitch;
    xall(evenRows) = xall(evenRows) + g.even_xOff;
    xall(oddRows) = xall(oddRows) + g.odd_xOff;

    yall = rowInd*g.vertPitch;
    
    figure('Name','shank view','Units','Normalized', 'Position', [0.2,0.1,0.3,0.8])
    for sI = 0:g.nShank-1
        cc = find(shankInd == sI);
        scatter( g.shankPitch*sI + xall, yall, 3, 'k', 'square' ); hold on;
        scatter( g.shankPitch*sI + xCoord(cc), yCoord(cc), 30, 'green', 'square', 'filled' ); hold on; 
    end
    xlim([min(xall)-5, (g.nShank-1)*g.shankPitch + max(xall)+5]);
    ylim([min(yall)-5, max(yall)+5]);
    hold off;
end % plotSaved


% =========================================================
% Build snsGeomMap from xy coordinates
%
function snsGeomStr = snsGeom(meta, shankInd, xCoord, yCoord, use, nShank, shankSep, shankWidth )
    % header
    % get probe part number; if absent, this is a 3A
    if isfield(meta,'imDatPrb_pn')
        pn = meta.imDatPrb_pn;
    else
        pn = '3A';
    end
    snsGeomStr = sprintf('~snsGeomMap=(%s,%d,%d,%d)', pn, nShank, shankSep, shankWidth);
    for i = 1:numel(shankInd)
        currEntry = sprintf('(%d:%g:%g:%d)', shankInd(i),xCoord(i),yCoord(i),use(i));
        snsGeomStr = sprintf('%s%s',snsGeomStr,currEntry);
    end
end % snsGeom


% =========================================================
% Parse snsGeomMap for XY coordinates
%
function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta)

    C = textscan(meta.snsGeomMap, '(%d:%d:%d:%d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
    shankInd = double(cell2mat(C(1)));
    xCoord = double(cell2mat(C(2)));
    yCoord = double(cell2mat(C(3)));
    connected = double(cell2mat(C(4)));

    % parse header for number of shanks
    geomStr = meta.snsGeomMap;
    headStr = extractBefore(geomStr,')(');
    headParts = split(headStr,',');
    nShank = str2double(headParts{2});
    shankWidth = str2double(headParts{4});
    shankPitch = str2double(headParts{3});
end % geomMapToGeom


% =========================================================
% Get XY coordinates from snsShankMap plus hard coded geom values
%
function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta)
    % get number of saved AP channels (some early metadata files have a
    % SYNC entry in the snsChanMap
    [nchan,~,~] = ChannelCountsIM(meta);

    C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
            'EndOfLine', ')', 'HeaderLines', 1 );
    shankInd = double(cell2mat(C(1)));
    colInd = double(cell2mat(C(2)));
    rowInd = double(cell2mat(C(3)));
    connected = double(cell2mat(C(4)));
    
    % trim these to the number of saved channels
    shankInd = shankInd(1:nchan);
    colInd = colInd(1:nchan);
    rowInd = rowInd(1:nchan);
    connected = connected(1:nchan);

    geom = getGeomParams(meta);

    oddRows = logical(mod(rowInd,2));
    evenRows = ~oddRows;
    xCoord = colInd*geom.horzPitch;
    xCoord(evenRows) = xCoord(evenRows) + geom.even_xOff ;
    xCoord(oddRows) = xCoord(oddRows) + geom.odd_xOff;
    yCoord = rowInd*geom.vertPitch;
    
    nShank = geom.nShank;
    shankWidth = geom.shankWidth;
    shankPitch = geom.shankPitch;
end % shankMapToGeom


% =========================================================
% Return geometry paramters for supported probe types
% These are used to calculate positions from metadata
% that includes only ~snsShankMap
%
function geom = getGeomParams(meta)
% create map
geomTypeMap = makeTypeMap();

% get probe part number; if absent, this is a 3A
if isfield(meta,'imDatPrb_pn')
    pn = meta.imDatPrb_pn;
else
    pn = '3A';
end

if geomTypeMap.isKey(pn)
    geomType = geomTypeMap(pn);
else
    fprintf('unsupported probe part number\n');
    return;
end

switch geomType
    case 'np1_stag_70um'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'nhp_lin_70um'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'nhp_stag_125um_med'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 87;
        geom.vertPitch = 20;
        geom.rowsPerShank = 1368;
        geom.elecPerShank = 2496;
    case 'nhp_stag_125um_long'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 87;
        geom.vertPitch = 20;
        geom.rowsPerShank = 2208;
        geom.elecPerShank = 4416;
    case 'nhp_lin_125um_med'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 103;
        geom.vertPitch = 20;
        geom.rowsPerShank = 1368;
        geom.elecPerShank = 2496;
    case 'nhp_lin_125um_long'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 103;
        geom.vertPitch = 20;
        geom.rowsPerShank = 2208;
        geom.elecPerShank = 4416;
    case 'uhd_8col_1bank'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 14;
        geom.odd_xOff = 14;
        geom.horzPitch = 6;
        geom.vertPitch = 6;
        geom.rowsPerShank = 48;
        geom.elecPerShank = 384;
   case 'uhd_8col_16bank'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 14;
        geom.odd_xOff = 14;
        geom.horzPitch = 6;
        geom.vertPitch = 6;
        geom.rowsPerShank = 768;
        geom.elecPerShank = 6144;
    case 'np2_ss'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 15;
        geom.rowsPerShank = 640;
        geom.elecPerShank = 1280;
     case 'np2_4s'
        geom.nShank = 4;
        geom.shankWidth = 70;
        geom.shankPitch = 250;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 15;
        geom.rowsPerShank = 640;
        geom.elecPerShank = 1280;
    case 'NP1120'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 6.75;
        geom.odd_xOff = 6.75;
        geom.horzPitch = 4.5;
        geom.vertPitch = 4.5;
        geom.rowsPerShank = 192;
        geom.elecPerShank = 384;
    case 'NP1121'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 6.25;
        geom.odd_xOff = 6.25;
        geom.horzPitch = 3;
        geom.vertPitch = 3;
        geom.rowsPerShank = 384;
        geom.elecPerShank = 384;
    case 'NP1122'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 12.5;
        geom.odd_xOff = 12.5;
        geom.horzPitch = 3;
        geom.vertPitch = 3;
        geom.rowsPerShank = 24;
        geom.elecPerShank = 384;
    case 'NP1123'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 10.25;
        geom.odd_xOff = 10.25;
        geom.horzPitch = 4.5;
        geom.vertPitch = 4.5;
        geom.rowsPerShank = 32;
        geom.elecPerShank = 384;
    case 'NP1300'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 48;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'NP1200'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 64;
        geom.elecPerShank = 128;
    case 'NXT3000'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 53;
        geom.odd_xOff = 53;
        geom.horzPitch = 0;
        geom.vertPitch = 15;
        geom.rowsPerShank = 128;
        geom.elecPerShank = 128;
    otherwise
        % shouldn't see this case
        fprintf('unsupported probe part number\n');
        return;
end
end % getGeomParams


% =========================================================
% Return geometry paramters for supported probe types
% Note that geom only contains enough info to calculate
% positions for the electrodes listed in snsShankMap
%
function M = makeTypeMap()
% many part numbers have the same geometry parameters ;
% make a map that pairs geometry type (value) with probe part number (key)
M = containers.Map('KeyType','char','ValueType','char');

M('3A') = 'np1_stag_70um';
M('PRB_1_4_0480_1') = 'np1_stag_70um';
M('PRB_1_4_0480_1_C') = 'np1_stag_70um';
M('NP1010') = 'np1_stag_70um'; 
M('NP1011') = 'np1_stag_70um';
M('NP1012') = 'np1_stag_70um';
M('NP1013') = 'np1_stag_70um';

M('NP1015') = 'nhp_lin_70um';
M('NP1015') = 'nhp_lin_70um';
M('NP1016') = 'nhp_lin_70um';
M('NP1017') = 'nhp_lin_70um';
   
M('NP1020') = 'nhp_stag_125um_med';
M('NP1021') = 'nhp_stag_125um_med';
M('NP1030') = 'nhp_stag_125um_long';
M('NP1031') = 'nhp_stag_125um_long';

M('NP1022') = 'nhp_lin_125um_med';
M('NP1032') = 'nhp_lin_125um_long';

M('NP1100') = 'uhd_8col_1bank';
M('NP1110') = 'uhd_8col_16bank';

M('PRB2_1_2_0640_0') = 'np2_ss';
M('PRB2_1_4_0480_1') = 'np2_ss';
M('NP2000') = 'np2_ss';
M('NP2003') = 'np2_ss';
M('NP2004') = 'np2_ss';

M('PRB2_4_2_0640_0') = 'np2_4s';
M('PRB2_4_4_0480_1') = 'np2_4s';
M('NP2010') = 'np2_4s';
M('NP2013') = 'np2_4s';
M('NP2014') = 'np2_4s';

M('NP1120') = 'NP1120';
M('NP1121') = 'NP1121';
M('NP1122') = 'NP1122';
M('NP1123') = 'NP1123';
M('NP1300') = 'NP1300';

M('NP1200') = 'NP1200';
M('NXT3000') = 'NXT3000';
end % makeTypeMap

% =========================================================
% Return full MUX table string to append to metadata
%
function muxTableStr = getMuxTable(meta)
% Read probe part number from meta
% Return full MUX table string to append to metadata
% As of 032923, there are 4 mux tables

np1='~muxTbl=(32,12)(0 1 24 25 48 49 72 73 96 97 120 121 144 145 168 169 192 193 216 217 240 241 264 265 288 289 312 313 336 337 360 361)(2 3 26 27 50 51 74 75 98 99 122 123 146 147 170 171 194 195 218 219 242 243 266 267 290 291 314 315 338 339 362 363)(4 5 28 29 52 53 76 77 100 101 124 125 148 149 172 173 196 197 220 221 244 245 268 269 292 293 316 317 340 341 364 365)(6 7 30 31 54 55 78 79 102 103 126 127 150 151 174 175 198 199 222 223 246 247 270 271 294 295 318 319 342 343 366 367)(8 9 32 33 56 57 80 81 104 105 128 129 152 153 176 177 200 201 224 225 248 249 272 273 296 297 320 321 344 345 368 369)(10 11 34 35 58 59 82 83 106 107 130 131 154 155 178 179 202 203 226 227 250 251 274 275 298 299 322 323 346 347 370 371)(12 13 36 37 60 61 84 85 108 109 132 133 156 157 180 181 204 205 228 229 252 253 276 277 300 301 324 325 348 349 372 373)(14 15 38 39 62 63 86 87 110 111 134 135 158 159 182 183 206 207 230 231 254 255 278 279 302 303 326 327 350 351 374 375)(16 17 40 41 64 65 88 89 112 113 136 137 160 161 184 185 208 209 232 233 256 257 280 281 304 305 328 329 352 353 376 377)(18 19 42 43 66 67 90 91 114 115 138 139 162 163 186 187 210 211 234 235 258 259 282 283 306 307 330 331 354 355 378 379)(20 21 44 45 68 69 92 93 116 117 140 141 164 165 188 189 212 213 236 237 260 261 284 285 308 309 332 333 356 357 380 381)(22 23 46 47 70 71 94 95 118 119 142 143 166 167 190 191 214 215 238 239 262 263 286 287 310 311 334 335 358 359 382 383)';
np2='~muxTbl=(24,16)(0 1 32 33 64 65 96 97 128 129 160 161 192 193 224 225 256 257 288 289 320 321 352 353)(2 3 34 35 66 67 98 99 130 131 162 163 194 195 226 227 258 259 290 291 322 323 354 355)(4 5 36 37 68 69 100 101 132 133 164 165 196 197 228 229 260 261 292 293 324 325 356 357)(6 7 38 39 70 71 102 103 134 135 166 167 198 199 230 231 262 263 294 295 326 327 358 359)(8 9 40 41 72 73 104 105 136 137 168 169 200 201 232 233 264 265 296 297 328 329 360 361)(10 11 42 43 74 75 106 107 138 139 170 171 202 203 234 235 266 267 298 299 330 331 362 363)(12 13 44 45 76 77 108 109 140 141 172 173 204 205 236 237 268 269 300 301 332 333 364 365)(14 15 46 47 78 79 110 111 142 143 174 175 206 207 238 239 270 271 302 303 334 335 366 367)(16 17 48 49 80 81 112 113 144 145 176 177 208 209 240 241 272 273 304 305 336 337 368 369)(18 19 50 51 82 83 114 115 146 147 178 179 210 211 242 243 274 275 306 307 338 339 370 371)(20 21 52 53 84 85 116 117 148 149 180 181 212 213 244 245 276 277 308 309 340 341 372 373)(22 23 54 55 86 87 118 119 150 151 182 183 214 215 246 247 278 279 310 311 342 343 374 375)(24 25 56 57 88 89 120 121 152 153 184 185 216 217 248 249 280 281 312 313 344 345 376 377)(26 27 58 59 90 91 122 123 154 155 186 187 218 219 250 251 282 283 314 315 346 347 378 379)(28 29 60 61 92 93 124 125 156 157 188 189 220 221 252 253 284 285 316 317 348 349 380 381)(30 31 62 63 94 95 126 127 158 159 190 191 222 223 254 255 286 287 318 319 350 351 382 383)';
np1100='~muxTbl=(32,12)(0 1 24 25 48 49 72 73 96 97 120 121 144 145 168 169 192 193 216 217 240 241 264 265 288 289 312 313 336 337 360 361)(2 3 26 27 50 51 74 75 98 99 122 123 146 147 170 171 194 195 218 219 242 243 266 267 290 291 314 315 338 339 362 363)(4 5 28 29 52 53 76 77 100 101 124 125 148 149 172 173 196 197 220 221 244 245 268 269 292 293 316 317 340 341 364 365)(6 7 30 31 54 55 78 79 102 103 126 127 150 151 174 175 198 199 222 223 246 247 270 271 294 295 318 319 342 343 366 367)(8 9 32 33 56 57 80 81 104 105 128 129 152 153 176 177 200 201 224 225 248 249 272 273 296 297 320 321 344 345 368 369)(10 11 34 35 58 59 82 83 106 107 130 131 154 155 178 179 202 203 226 227 250 251 274 275 298 299 322 323 346 347 370 371)(12 13 36 37 60 61 84 85 108 109 132 133 156 157 180 181 204 205 228 229 252 253 276 277 300 301 324 325 348 349 372 373)(14 15 38 39 62 63 86 87 110 111 134 135 158 159 182 183 206 207 230 231 254 255 278 279 302 303 326 327 350 351 374 375)(16 17 40 41 64 65 88 89 112 113 136 137 160 161 184 185 208 209 232 233 256 257 280 281 304 305 328 329 352 353 376 377)(18 19 42 43 66 67 90 91 114 115 138 139 162 163 186 187 210 211 234 235 258 259 282 283 306 307 330 331 354 355 378 379)(20 21 44 45 68 69 92 93 116 117 140 141 164 165 188 189 212 213 236 237 260 261 284 285 308 309 332 333 356 357 380 381)(22 23 46 47 70 71 94 95 118 119 142 143 166 167 190 191 214 215 238 239 262 263 286 287 310 311 334 335 358 359 382 383)';
np128ch='~muxTbl=(12,12)(84 11 85 5 74 10 56 112 46 121 39 127)(100 26 110 33 69 24 63 109 45 93 25 99)(87 0 82 6 71 15 53 117 43 122 42 116)(102 28 81 34 70 18 60 103 17 94 27 101)(73 1 86 7 68 16 50 106 40 123 128 128)(105 29 75 35 67 12 54 89 20 95 128 128)(76 2 83 8 65 13 47 118 49 124 128 128)(108 30 78 36 64 14 51 90 23 96 128 128)(79 3 80 9 62 114 44 119 52 125 128 128)(104 31 72 37 61 113 57 91 19 97 128 128)(88 4 77 21 59 111 41 120 55 126 128 128)(107 32 66 38 58 115 48 92 22 98 128 128)';

M = containers.Map('KeyType','char','ValueType','char');

M('3A') = np1;
M('PRB_1_4_0480_1') = np1;
M('PRB_1_4_0480_1_C') = np1;
M('NP1010') = np1; 
M('NP1011') = np1;
M('NP1012') = np1;
M('NP1013') = np1;

M('NP1015') = np1;
M('NP1015') = np1;
M('NP1016') = np1;
M('NP1017') = np1;
   
M('NP1020') = np1;
M('NP1021') = np1;
M('NP1030') = np1;
M('NP1031') = np1;

M('NP1022') = np1;
M('NP1032') = np1;

M('NP1100') = np1;
M('NP1110') = np1100;

M('PRB2_1_4_0480_1') = np2;
M('PRB2_1_2_0640_0') = np2;
M('NP2000') = np2;
M('NP2003') = np2;
M('NP2004') = np2;

M('PRB2_4_2_0640_0') = np2;
M('PRB2_4_4_0480_1') = np2;
M('NP2010') = np2;
M('NP2013') = np2;
M('NP2014') = np2;

M('NP1120') = np1;
M('NP1121') = np1;
M('NP1122') = np1;
M('NP1123') = np1;
M('NP1300') = np1;

M('NP1200') = np128ch;
M('NXT3000') = np128ch;

% get probe part number; if absent, this is a 3A
if isfield(meta,'imDatPrb_pn')
    pn = meta.imDatPrb_pn;
else
    pn = '3A';
end

if M.isKey(pn)
    muxTableStr = M(pn);
else
    muxTableStr = '';
    fprintf('unsupported probe part number\n');
    return;
end   

end


% =========================================================
% Parse imro table to extract 'new' metadata items
% (SpikeGLX 032623-phase 30 and later)
% Returns strings for:
%     imChan0apGain
%     imChan0lfGain
%     imAnyChanFullBand
% 
%
function [imChan0apGainStr, imChan0lfGainStr, imAnyChanFullBandStr] = imroMetaItems(meta)
   
% read header to get what type of imro this is
headEntry = extractAfter(extractBefore(meta.imroTbl, ')'), '(');
headList = split(headEntry, ',');
currType = headList{1};

if str2num(currType) > 50000
    % this is a 3A probe
    currType = '3A';  % in matlab, easier to parse imroTbl specifically for 3A
end

switch currType

    case {'24', '21'}
        imChan0apGainStr = 'imChan0apGain=80';
        imChan0lfGainStr = 'imChan0lfGain=80';
        imAnyChanFullBandStr = 'imAnyChanFullBand=true';
    
    case '0'
        C = textscan(meta.imroTbl, '(%d %d %d %d %d %d', 'EndOfLine', ')', 'HeaderLines', 1 );
        % parse first entry of the table for ap and lf gain
        imChan0apGainStr = sprintf('imChan0apGain=%d', C{4}(1));
        imChan0lfGainStr = sprintf('imChan0lfGain=%d', C{5}(1));
        % get array of APFilt settings
        APFilt = double(cell2mat(C(6))); % array of all APFilt settings
        if sum(APFilt==0) > 0     
            imAnyChanFullBandStr = 'imAnyChanFullBand=true';
        else
            imAnyChanFullBandStr = 'imAnyChanFullBand=false';
        end
     
    case '3A'
        C = textscan(meta.imroTbl, '(%d %d %d %d %d', 'EndOfLine', ')', 'HeaderLines', 1 );
        % parse first entry of the table for ap and lf gain
        imChan0apGainStr = sprintf('imChan0apGain=%d', C{4}(1));
        imChan0lfGainStr = sprintf('imChan0lfGain=%d', C{5}(1));
        % no APFilt setting available
        imAnyChanFullBandStr = 'imAnyChanFullBand=false';
    
    case '1110'
        % ap and lf gain and filter option are in the header 
        imChan0apGainStr = sprintf('imChan0apGain=%s', headList{4});
        imChan0lfGainStr = sprintf('imChan0lfGain=%s', headList{5});
        if strcmp(headList{6}, '0')
            imAnyChanFullBandStr = 'imAnyChanFullBand=true';
        else
            imAnyChanFullBandStr = 'imAnyChanFullBand=false';
        end
end
   
end