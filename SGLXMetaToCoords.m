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
%     3 creates a new metadata file with an snsGeomMap
%
% File is saved in current working directory 
%
% Jennifer Colonell, Janelia Research Campus
 
function SGLXMetaToCoords()

% Output selection:
outType = 1;

% Plot shank view
bPlot = 1;

% Ask user for metadata file
[metaName,path] = uigetfile('*.meta', 'Select Metadata File');

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

    case 3 % new metadata file that includes the snsGeomMap field
        newName = [fname,'_plusGeom.meta'];
        copyfile(fullfile(path,metaName), fullfile(path,newName));
        fid = fopen(fullfile(path,newName), 'a');
        snsGeomStr = snsGeom(meta, shankInd, xCoord, yCoord, connected, nShank, shankPitch, shankWidth );
        fprintf(fid,'%s\n',snsGeomStr);
        fclose(fid);
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
    snsGeomStr = sprintf('~snsGeomMap=(%s,%d,%d,%d)', meta.imDatPrb_pn, nShank, shankSep, shankWidth);
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


