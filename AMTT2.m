clear
clc
close all

%% instructions
% The input file containing aspect data must be in TIFF or GeoTIFF format.
% The input file containing slope data must be in TIFF or GeoTIFF format.
% The input file containing orientation of discontinuity sets must contain, in this exact order: [dip, dip direction, friction angle]. It must be in xls or xlsx format.
% Output files are saved in DIR_output directory. They consist of TIFF or GeoTIFF files (same format as the one of the input slope and aspect files) that can be directly input in Gis environment.

%% names of output files
% Deterministic mode
fileout2D_D = 'D_2D_';% prefix of files containing results of planar (2D) sliding (it will be followed by the id of the considered set and 'sum' in the case of sum of all possible planar sliding kinematisms)
fileout3D_D = 'D_3D_';% prefix of files containing results of wedge (3D) sliding (it will be followed by the id of the considered combination of sets and 'sum' in the case of sum of all possible wedge sliding kinematisms)
fileoutDT_D = 'D_DT_';% prefix of files containing results of direct toppling (it will be followed by the id of the considered combination of sets and 'sum' in the case of sum of all possible direct toppling kinematisms)
fileoutFT_D = 'D_FT_';% prefix of files containing results of flexural toppling (it will be followed by the id of the considered set and 'sum' in the case of sum of all possible flexural toppling kinematisms)
fileoutfinal_D = 'D_RESULTS';% file containing results of all the analized phenomena
fileoutsource_D = 'D_SOURCES';% file containing source areas

% Raw Data mode (Custom)
fileout2D_C = 'C_2D_';% prefix of files containing results of planar (2D) sliding (followed by number of the considered set or 'sum' in case of global result)
fileout3D_C = 'C_3D_';% prefix of files containing results of wedge (3D) sliding (followed by number of the two considered sets or 'sum' in case of global result)
fileoutDT_C = 'C_DT_';% prefix of files containing results of direct toppling (followed by number of the two considered sets or 'sum' in case of global result)
fileoutFT_C = 'C_FT_';% prefix of files containing results of flexural toppling (followed by number of the considered set or 'sum' in case of global result)
fileoutfinal_C = 'C_RESULTS';% file containing results of all the analized kinematisms
fileoutsource_C = 'C_SOURCES';% file containing source areas

%% DO NOT CHANGE THE CODE FROM HERE ON OUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT DIRECTORY %%%%%%%%%%%%%%
[aspectfile,path_aspect] = uigetfile({'*.tiff';'*.tif'},'Select the input aspect file');
[slopefile,path_slope] = uigetfile({'*.tiff';'*.tif'},'Select the input slope file');
[orientationfile,path_orientation] = uigetfile({'*.xlsx';'*.xls'},'Select the input orientation file');
aspect = imread([path_aspect aspectfile]);
info_aspect = geotiffinfo([path_aspect aspectfile]);
slope = imread([path_slope slopefile]);
info_slope = geotiffinfo([path_slope slopefile]);
Orientazione = readtable([path_orientation orientationfile]);
Orientazione = table2array(Orientazione);

%%%% OUTPUT DIRECTORY %%%%%%%%%%%%%%
DIR_output = uigetdir('C:\', 'Choose the output directory');

%%%% USER INTERFACE %%%%%%%%%%%%%%%%
fig = uifigure('Name','Choose among methods for performing kinematic analysis','Position',[400 400 750 300]);

% botton group "kinematic analysis mode"
bg = uibuttongroup(fig,'Position',[20 100 120 135]);
tb1 = uitogglebutton(bg,'Position',[10 75 100 22],'Text','Deterministic');
tb3 = uitogglebutton(bg,'Position',[10 25 100 22],'Text','Raw Data');
c1 = uicontrol(fig,'Position', [30 70 100 22],'String','Confirm method','Callback','uiresume(fig)');
uiwait(fig)

% botton "Start calculation"--- START OF AUTOMATIC OPERATIONS
c2 = uicontrol(fig,'Position', [30 20 100 22],'String','Start calculation','Callback','uiresume(fig)');
uiwait(fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tb1.Value == true
    mode = 1;% Deterministic

    %% discontinuity sets
    phi = min(Orientazione(:,3));% minimum value of phi
    nset = size(Orientazione,1);

    %% operations
    aspect(aspect==-9999)= NaN;
    slope(slope==-9999)= NaN;

    VER2Dsum = zeros(size(aspect,1),size(aspect,2));
    VER3Dsum = zeros(size(aspect,1),size(aspect,2));
    VERFTsum = zeros(size(aspect,1),size(aspect,2));
    VERDTsum = zeros(size(aspect,1),size(aspect,2));
    VERsum = zeros(size(aspect,1),size(aspect,2));

    %% MARKLAND'S TEST FOR PLANAR (2D) SLIDING

    for i = 1:nset

        % riclassification of aspect
        M_dipdir = Orientazione(i,2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - M_dipdir);
        Tmax = 40;
        Tmin = 0;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmax)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmin)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Tphi = Orientazione(i,3);
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<Tphi)= NaN;
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        Verifica2D = Condizione_orientazione_reclass & Condizione_pendenza;

        Verifica2D = double(Verifica2D);
        Verifica2Dnonriclass = Verifica2D;

        % write geotiff in output
        filenameout = [fileout2D_D num2str(i)];
        geotiffwrite([DIR_output '\' filenameout],Verifica2D,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VER2D(:,:,i) = Verifica2D;
        VER2Dsum = VER2Dsum + Verifica2Dnonriclass;
    end

    filenameout = [fileout2D_D 'sum'];
    geotiffwrite([DIR_output '\' filenameout],VER2Dsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);


    %% MARKLAND'S TEST FOR WEDGE (3D) SLIDING

    % identification of all the possible combinations of sets
    n_col = 2;
    t = 1:nset;
    y = repmat(t,1,nset^(n_col-1))';
    x = repmat(t,nset,nset^(n_col-2));
    xr = reshape(x,[size(y) 1]);
    C = [y xr];
    comb = sort(C,2);
    comb = unique(comb, 'rows');
    sy = [1:size(comb(:,1))]';
    ind = comb(:,1)-comb(:,2);
    ind(:,2) = sy;
    ind(ind==0,:) = [];
    Comb_set = comb(ind(:,2),:);

    for j = 1:size(Comb_set,1)
        dip1 = Orientazione(Comb_set(j,1),1);
        dipdir1 = Orientazione(Comb_set(j,1),2);
        dip2 = Orientazione(Comb_set(j,2),1);
        dipdir2 = Orientazione(Comb_set(j,2),2);
        n1 = [sin(dip1*pi/180)*sin(dipdir1*pi/180) sin(dip1*pi/180)*cos(dipdir1*pi/180) cos(dip1*pi/180)];
        n2 = [sin(dip2*pi/180)*sin(dipdir2*pi/180) sin(dip2*pi/180)*cos(dipdir2*pi/180) cos(dip2*pi/180)];
        l1 = norm(n1);
        l2 = norm(n2);
        pv = cross(n1/l1,n2/l2);

        if pv(3) < 0
            pv = -pv;
        end
        dip = acos(abs(pv(3))/norm(pv,'fro'))*180/pi;

        angle = acos(pv(2)/norm(pv(1:2),'fro'))*180/pi;
        if pv(3) >= 0
            if pv(1) >= 0
                dipDir = angle;
            else
                dipDir = 360-angle;
            end
        else
            if pv(1) >= 0
                dipDir = 180 + angle;
            else
                dipDir = 180 - angle;
            end
        end
        intersezione = [dip dipDir];

        % riclassification of aspect
        I_dipdir = intersezione(2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - I_dipdir);
        Tmaxi = 180;
        Tmini = 0;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmaxi)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmini)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<intersezione(1))= NaN;% questa riga è modificata rispetto all'originale, perchè mi sembrava sbagliato il >: verificare!
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        Verifica3D = Condizione_orientazione_reclass & Condizione_pendenza;

        Verifica3D = double(Verifica3D);
        Verifica3Dnonriclass = Verifica3D;

        % write geotiff in output
        filenameout = [fileout3D_D num2str(Comb_set(j,1)) '-' num2str(Comb_set(j,2))];
        geotiffwrite([DIR_output '\' filenameout],Verifica3D,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VER3D(:,:,j) = Verifica3D;
        VER3Dsum = VER3Dsum + Verifica3Dnonriclass;
    end

    filenameout = [fileout3D_D 'sum'];
    geotiffwrite([DIR_output '\' filenameout],VER3Dsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);


    %% MARKLAND'S TEST FOR DIRECT TOPPLING

    % identification of all the possible combinations of sets
    n_col = 2;
    t = 1:nset;
    y = repmat(t,1,nset^(n_col-1))';
    x = repmat(t,nset,nset^(n_col-2));
    xr = reshape(x,[size(y) 1]);
    C = [y xr];
    comb = sort(C,2);
    comb = unique(comb, 'rows');
    sy = [1:size(comb(:,1))]';
    ind = comb(:,1)-comb(:,2);
    ind(:,2) = sy;
    ind(ind==0,:) = [];
    Comb_set = comb(ind(:,2),:);

    for j = 1:size(Comb_set,1)
        dip1 = Orientazione(Comb_set(j,1),1);
        dipdir1 = Orientazione(Comb_set(j,1),2);
        dip2 = Orientazione(Comb_set(j,2),1);
        dipdir2 = Orientazione(Comb_set(j,2),2);
        n1 = [sin(dip1*pi/180)*sin(dipdir1*pi/180) sin(dip1*pi/180)*cos(dipdir1*pi/180) cos(dip1*pi/180)];
        n2 = [sin(dip2*pi/180)*sin(dipdir2*pi/180) sin(dip2*pi/180)*cos(dipdir2*pi/180) cos(dip2*pi/180)];
        l1 = norm(n1);
        l2 = norm(n2);
        pv = cross(n1/l1,n2/l2);

        if pv(3) < 0
            pv = -pv;
        end
        dip = acos(abs(pv(3))/norm(pv,'fro'))*180/pi;

        angle = acos(pv(2)/norm(pv(1:2),'fro'))*180/pi;
        if pv(3) >= 0
            if pv(1) >= 0
                dipDir = angle;
            else
                dipDir = 360-angle;
            end
        else
            if pv(1) >= 0
                dipDir = 180 + angle;
            else
                dipDir = 180 - angle;
            end
        end
        intersezione = [dip dipDir];

        % riclassification of aspect
        I_dipdir = intersezione(2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - I_dipdir);
        Tmaxi = 40;
        Tmini = 0;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmaxi)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmini)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<intersezione(1))= NaN;% stessa modifica fatta per 3D
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        VerificaDT = Condizione_orientazione_reclass & Condizione_pendenza;

        VerificaDT = double(VerificaDT);
        VerificaDTnonriclass = VerificaDT;

        % write geotiff in output
        filenameout = [fileoutDT_D num2str(Comb_set(j,1)) '-' num2str(Comb_set(j,2))];
        geotiffwrite([DIR_output '\' filenameout],VerificaDT,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VERDT(:,:,j) = VerificaDT;
        VERDTsum = VERDTsum + VerificaDTnonriclass;
    end

    filenameout = [fileoutDT_D 'sum'];
    geotiffwrite([DIR_output '\' filenameout],VERDTsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% MARKLAND'S TEST FOR FLEXURAL TOPPLING

    for i = 1:nset

        % riclassification of aspect
        M_dipdir = Orientazione(i,2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - M_dipdir);
        Tmax = 200;
        Tmin = 160;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmax)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmin)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Ttop = (90 - Orientazione(i,1))+ Orientazione(i,3);
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<Ttop)= NaN;
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        VerificaFT = Condizione_orientazione_reclass & Condizione_pendenza;

        VerificaFT = double(VerificaFT);
        VerificaFTnonriclass = VerificaFT;

        % write geotiff in output
        filenameout = [fileoutFT_D num2str(i)];
        geotiffwrite([DIR_output '\' filenameout],VerificaFT,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VERFT(:,:,i) = VerificaFT;
        VERFTsum = VERFTsum + VerificaFTnonriclass;
    end

    filenameout = [fileoutFT_D 'sum'];
    geotiffwrite([DIR_output '\' filenameout],VERFTsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% SUM OF THE RESULTS OF ALL THE POSSIBLE KINEMATISMS
    VERsum = VER2Dsum + VER3Dsum + VERDTsum + VERFTsum;
    geotiffwrite([DIR_output '\' fileoutfinal_D],VERsum,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% CREATION OF THE RASTER FILE CONTAINING THE SOURCE POINTS
    VERsource = VERsum;
    VERsource(VERsource>1)= 1;
    geotiffwrite([DIR_output '\' fileoutsource_D],VERsource,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% end of the process
    sprintf('Calculation ended. Enjoy!')
    close(fig)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif tb3.Value == true
    mode = 3;% Custom (Raw Data)

    %% discontinuity planes
    phi = min(Orientazione(:,3));% minimum value of phi
    nplanes = size(Orientazione,1);

    %% operations
    aspect(aspect==-9999)= NaN;
    slope(slope==-9999)= NaN;

    VER2Dsum = zeros(size(aspect,1),size(aspect,2));
    VER3Dsum = zeros(size(aspect,1),size(aspect,2));
    VERFTsum = zeros(size(aspect,1),size(aspect,2));
    VERDTsum = zeros(size(aspect,1),size(aspect,2));
    VERsum = zeros(size(aspect,1),size(aspect,2));

    %% MARKLAND'S TEST FOR PLANAR (2D) SLIDING

    for i = 1:nplanes

        % riclassification of aspect
        M_dipdir = Orientazione(i,2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - M_dipdir);
        Tmax = 40;
        Tmin = 0;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmax)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmin)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Tphi = Orientazione(i,3);
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<Tphi)= NaN;
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        Verifica2D = Condizione_orientazione_reclass & Condizione_pendenza;

        Verifica2D = double(Verifica2D);
        Verifica2Dnonriclass = Verifica2D;

        % write geotiff in output
        filenameout = [fileout2D_C num2str(i)];
        DIR_output2D = [DIR_output '\' fileout2D_C(1:end-1)];
        if i == 1
            mkdir(DIR_output2D);
        end
        geotiffwrite([DIR_output2D '\' filenameout],Verifica2D,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VER2D(:,:,i) = Verifica2D;
        VER2Dsum = VER2Dsum + Verifica2Dnonriclass;
    end


    VER2Dtest = VER2Dsum/nplanes;
    filenameout = [fileout2D_C 'sum'];
    geotiffwrite([DIR_output '\' filenameout],VER2Dtest,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);


    %% MARKLAND'S TEST FOR WEDGE (3D) SLIDING

    % identification of all the possible combinations of sets
    n_col = 2;
    t = 1:nplanes;
    y = repmat(t,1,nplanes^(n_col-1))';
    x = repmat(t,nplanes,nplanes^(n_col-2));
    xr = reshape(x,[size(y) 1]);
    C = [y xr];
    comb = sort(C,2);
    comb = unique(comb, 'rows');
    sy = [1:size(comb(:,1))]';
    ind = comb(:,1)-comb(:,2);
    ind(:,2) = sy;
    ind(ind==0,:) = [];
    Comb_planes = comb(ind(:,2),:);
    ncomb3D = size(Comb_planes,1);

    for j = 1:ncomb3D
        dip1 = Orientazione(Comb_planes(j,1),1);
        dipdir1 = Orientazione(Comb_planes(j,1),2);
        dip2 = Orientazione(Comb_planes(j,2),1);
        dipdir2 = Orientazione(Comb_planes(j,2),2);
        n1 = [sin(dip1*pi/180)*sin(dipdir1*pi/180) sin(dip1*pi/180)*cos(dipdir1*pi/180) cos(dip1*pi/180)];
        n2 = [sin(dip2*pi/180)*sin(dipdir2*pi/180) sin(dip2*pi/180)*cos(dipdir2*pi/180) cos(dip2*pi/180)];
        l1 = norm(n1);
        l2 = norm(n2);
        pv = cross(n1/l1,n2/l2);

        if pv(3) < 0
            pv = -pv;
        end
        dip = acos(abs(pv(3))/norm(pv,'fro'))*180/pi;

        angle = acos(pv(2)/norm(pv(1:2),'fro'))*180/pi;
        if pv(3) >= 0
            if pv(1) >= 0
                dipDir = angle;
            else
                dipDir = 360-angle;
            end
        else
            if pv(1) >= 0
                dipDir = 180 + angle;
            else
                dipDir = 180 - angle;
            end
        end
        intersezione = [dip dipDir];

        % riclassification of aspect
        I_dipdir = intersezione(2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - I_dipdir);
        Tmaxi = 180;
        Tmini = 0;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmaxi)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmini)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<intersezione(1))= NaN;% questa riga è modificata rispetto all'originale, perchè mi sembrava sbagliato il >: verificare!
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        Verifica3D = Condizione_orientazione_reclass & Condizione_pendenza;

        Verifica3D = double(Verifica3D);
        Verifica3Dnonriclass = Verifica3D;

        % write geotiff in output
        filenameout = [fileout3D_C num2str(Comb_planes(j,1)) '-' num2str(Comb_planes(j,2))];
        DIR_output3D = [DIR_output '\' fileout3D_C(1:end-1)];
        if j == 1
            mkdir(DIR_output3D);
        end
        geotiffwrite([DIR_output3D '\' filenameout],Verifica3D,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VER3D(:,:,j) = Verifica3D;
        VER3Dsum = VER3Dsum + Verifica3Dnonriclass;
    end

    VER3Dtest = VER3Dsum/ncomb3D;
    filenameout = [fileout3D_C 'sum'];
    geotiffwrite([DIR_output '\' filenameout],VER3Dtest,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);


    %% TEST DI MARKLAND PER IL RIBALTAMENTO DIRETTO

    % identification of all the possible combinations of sets
    n_col = 2;
    t = 1:nplanes;
    y = repmat(t,1,nplanes^(n_col-1))';
    x = repmat(t,nplanes,nplanes^(n_col-2));
    xr = reshape(x,[size(y) 1]);
    C = [y xr];
    comb = sort(C,2);
    comb = unique(comb, 'rows');
    sy = [1:size(comb(:,1))]';
    ind = comb(:,1)-comb(:,2);
    ind(:,2) = sy;
    ind(ind==0,:) = [];
    Comb_planes = comb(ind(:,2),:);
    ncombDT = size(Comb_planes,1);

    for j = 1:ncombDT
        dip1 = Orientazione(Comb_planes(j,1),1);
        dipdir1 = Orientazione(Comb_planes(j,1),2);
        dip2 = Orientazione(Comb_planes(j,2),1);
        dipdir2 = Orientazione(Comb_planes(j,2),2);
        n1 = [sin(dip1*pi/180)*sin(dipdir1*pi/180) sin(dip1*pi/180)*cos(dipdir1*pi/180) cos(dip1*pi/180)];
        n2 = [sin(dip2*pi/180)*sin(dipdir2*pi/180) sin(dip2*pi/180)*cos(dipdir2*pi/180) cos(dip2*pi/180)];
        l1 = norm(n1);
        l2 = norm(n2);
        pv = cross(n1/l1,n2/l2);

        if pv(3) < 0
            pv = -pv;
        end
        dip = acos(abs(pv(3))/norm(pv,'fro'))*180/pi;

        angle = acos(pv(2)/norm(pv(1:2),'fro'))*180/pi;
        if pv(3) >= 0
            if pv(1) >= 0
                dipDir = angle;
            else
                dipDir = 360-angle;
            end
        else
            if pv(1) >= 0
                dipDir = 180 + angle;
            else
                dipDir = 180 - angle;
            end
        end
        intersezione = [dip dipDir];

        % riclassification of aspect
        I_dipdir = intersezione(2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - I_dipdir);
        Tmaxi = 40;
        Tmini = 0;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmaxi)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmini)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<intersezione(1))= NaN;% stessa modifica fatta per 3D
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        VerificaDT = Condizione_orientazione_reclass & Condizione_pendenza;

        VerificaDT = double(VerificaDT);
        VerificaDTnonriclass = VerificaDT;

        % write geotiff in output
        filenameout = [fileoutDT_C num2str(Comb_planes(j,1)) '-' num2str(Comb_planes(j,2))];
        DIR_outputDT = [DIR_output '\' fileoutDT_C(1:end-1)];
        if j == 1
            mkdir(DIR_outputDT);
        end
        geotiffwrite([DIR_outputDT '\' filenameout],VerificaDT,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VERDT(:,:,j) = VerificaDT;
        VERDTsum = VERDTsum + VerificaDTnonriclass;
    end

    VERDTtest = VERDTsum/ncombDT;
    filenameout = [fileoutDT_C 'sum'];
    geotiffwrite([DIR_output '\' filenameout],VERDTtest,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% TEST DI MARKLAND PER IL RIBALTAMENTO FLESSIONALE

    for i = 1:nplanes

        % riclassification of aspect
        M_dipdir = Orientazione(i,2)*ones(size(aspect,1),size(aspect,2));
        Condizione_orientazione = abs(aspect - M_dipdir);
        Tmax = 200;
        Tmin = 160;

        Condizione_orientazione_reclass = Condizione_orientazione;

        Condizione_orientazione_reclass(Condizione_orientazione>Tmax)= NaN;
        Condizione_orientazione_reclass(Condizione_orientazione<Tmin)= NaN;

        iNan = isnan(Condizione_orientazione_reclass);
        Condizione_orientazione_reclass(iNan==0) = 1;
        Condizione_orientazione_reclass(iNan==1) = 0;
        Condizione_orientazione_reclass = double(Condizione_orientazione_reclass);

        % riclassification of slope
        Ttop = (90 - Orientazione(i,1))+ Orientazione(i,3);
        Condizione_pendenza = slope;
        Condizione_pendenza(Condizione_pendenza<Ttop)= NaN;
        iNan = isnan(Condizione_pendenza);
        Condizione_pendenza(iNan==0) = 1;
        Condizione_pendenza(iNan==1) = 0;
        Condizione_pendenza = double(Condizione_pendenza);

        % kinematic test
        VerificaFT = Condizione_orientazione_reclass & Condizione_pendenza;

        VerificaFT = double(VerificaFT);
        VerificaFTnonriclass = VerificaFT;

        % write geotiff in output
        filenameout = [fileoutFT_C num2str(i)];
        DIR_outputFT = [DIR_output '\' fileoutFT_C(1:end-1)];
        if i == 1
            mkdir(DIR_outputFT);
        end
        geotiffwrite([DIR_outputFT '\' filenameout],VerificaFT,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

        VERFT(:,:,i) = VerificaFT;
        VERFTsum = VERFTsum + VerificaFTnonriclass;
    end

    VERFTtest = VERFTsum/nplanes;
    filenameout = [fileoutFT_C 'sum'];
    geotiffwrite([DIR_output '\' filenameout],VERFTtest,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% SUM OF THE RESULTS OF ALL THE POSSIBLE KINEMATISMS
    Ntot_test = nplanes + ncomb3D + ncombDT + nplanes;
    VERtest = 0.25*VER2Dtest + 0.25*VER3Dtest + 0.25*VERDTtest + 0.25*VERFTtest;
    geotiffwrite([DIR_output '\' fileoutfinal_C],VERtest,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% CREATION OF THE RASTER FILE CONTAINING THE SOURCE POINTS
    VERsource = VERtest;
    VERsource(VERsource>0)= 1;
    geotiffwrite([DIR_output '\' fileoutsource_C],VERsource,info_aspect.RefMatrix,'CoordRefSysCode',info_aspect.GeoTIFFCodes.PCS);

    %% end of the process
    sprintf('Calculation ended. Enjoy!')
    close(fig)

end
