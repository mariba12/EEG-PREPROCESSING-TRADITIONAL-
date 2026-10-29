%% 01_CREAR_PSD_LONG_EEG_CREA_SIN_OC.m
% EEG_CREA: unifica las PSD ORIGINALES (sin limpiar ni ponderar).
% Salida CSV, porque la tabla completa supera 1.048.576 filas de Excel.
% NO calcula specparam, NO colapsa ROIs y NO elimina atipicos.
% Electrodo numerico 1-19 segun correspondencia del estudio.
%
% MOMENTOS NUMERICOS:
%   OA=1; tSMS_1..tSMS_15 = 2..16; POST_1..POST_5 = 17..21.
%   OC esta EXCLUIDO y no recibe codigo.
%
% Cada libro Sxx.xlsx: hojas Fp1 ... Oz; columnas Frequencies, OA, OC,
% tSMS_1 ... tSMS_15, POST_1 ... POST_5.
% El script NO importa ni exporta la columna OC.
%
% Autor: EEG_CREA | Auditoria previa a specparam

clearvars; clc;

%% RUTAS
PSD_DIR = 'C:\Users\evadi\OneDrive\Escritorio\FENNSI\CREATIVIDAD_DEFINITIVOS_EEG\PSD_archivos_brutos_por_sujeto';
OUTPUT_DIR = 'C:\Users\evadi\OneDrive\Escritorio\FENNSI\CREATIVIDAD_DEFINITIVOS_EEG\Database_Resultados\PSD';
AUDIT_DIR = OUTPUT_DIR;
if ~isfolder(AUDIT_DIR), mkdir(AUDIT_DIR); end

% Nuevos nombres para NO sobrescribir el CSV original con OC.
OUT_LONG = fullfile(AUDIT_DIR,'PSD_LONG_2MIN_SIN_OC.csv');
OUT_QC   = fullfile(AUDIT_DIR,'QC_IMPORTACION_PSD_LONG_SIN_OC.csv');
OUT_MAP  = fullfile(AUDIT_DIR,'EQUIVALENCIA_ELECTRODOS.csv');
OUT_INCID = fullfile(AUDIT_DIR,'INCIDENCIAS_IMPORTACION_SIN_OC.csv');

%% PARAMETROS (sin umbrales de rechazo)
FMIN = 0.5;
FMAX = 40;
ELECTRODOS = {'Fp1','Fp2','F4','F7','F8','Fz', ...
              'C3','C4','Cz','T3','T4','T5','T6', ...
              'P3','P4','Pz','O1','O2','Oz'};

% Solo OA, estimulacion y POST. OC se omite por completo.
MOMENTOS = [{'OA'}, ...
            arrayfun(@(k)sprintf('tSMS_%d',k),1:15,'UniformOutput',false), ...
            arrayfun(@(k)sprintf('POST_%d',k),1:5,'UniformOutput',false)];
NOM_MOMENTOS = string(MOMENTOS);
CODIGOS_MOMENTOS = (1:21)';
N_MOM = numel(MOMENTOS);

INICIOS = [NaN,0:2:28,0:2:8];
FINES   = [NaN,2:2:30,2:2:10];
FASES = ["PRE_OA",repmat("tSMS",1,15),repmat("POST",1,5)];

assert(N_MOM==21 && numel(CODIGOS_MOMENTOS)==N_MOM);
assert(numel(INICIOS)==N_MOM && numel(FINES)==N_MOM && numel(FASES)==N_MOM);
assert(isequal(CODIGOS_MOMENTOS,(1:21)'));

% El indice de ELECTRODOS coincide con el codigo 1..19.
mapa = table((1:numel(ELECTRODOS))',string(ELECTRODOS(:)), ...
    'VariableNames',{'Electrodo','NombreElectrodo'});
writetable(mapa,OUT_MAP);

%% BUSCAR ARCHIVOS: SOLO Sxx.xlsx, NO LIBROS DE LONGITUD
dd = dir(fullfile(PSD_DIR,'S*.xlsx'));
esSujeto = ~cellfun(@isempty,regexp({dd.name},'^S\d+\.xlsx$','once'));
dd = dd(esSujeto);
if isempty(dd)
    error('No encuentro archivos Sxx.xlsx directamente en: %s',PSD_DIR);
end
[~,ix] = sort({dd.name});
dd = dd(ix);
fprintf('Libros PSD encontrados: %d\n',numel(dd));

% Borrar UNICAMENTE las salidas SIN_OC de una ejecucion previa.
if isfile(OUT_LONG), delete(OUT_LONG); end
if isfile(OUT_QC), delete(OUT_QC); end
if isfile(OUT_INCID), delete(OUT_INCID); end

primera = true;
incidencias = strings(0,1);
sujetoQC = strings(0,1);
canalQC = zeros(0,1);
momentoQC = zeros(0,1);
nFrecuenciasQC = zeros(0,1);
nValidosQC = zeros(0,1);
nInvalidosQC = zeros(0,1);
nAusentesQC = zeros(0,1);

totalFilas = 0;
for iS = 1:numel(dd)
    archivo = fullfile(dd(iS).folder,dd(iS).name);
    sujeto = string(erase(dd(iS).name,'.xlsx'));
    fprintf('Procesando %s (%d/%d)\n',sujeto,iS,numel(dd));

    tablasSujeto = cell(numel(ELECTRODOS),1);
    for iE = 1:numel(ELECTRODOS)
        canal = ELECTRODOS{iE};
        try
            T = readtable(archivo,'Sheet',canal,'VariableNamingRule','preserve');
        catch ME
            incidencias(end+1,1) = sujeto + " | " + canal + ...
                " | No se pudo leer hoja: " + string(ME.message);
            continue;
        end

        nombres = string(T.Properties.VariableNames);
        ixFreq = find(strcmpi(strtrim(nombres),"Frequencies"),1);
        if isempty(ixFreq)
            incidencias(end+1,1) = sujeto + " | " + canal + " | Falta columna Frequencies";
            continue;
        end
        freq = convertirNumerico(T{:,ixFreq});
        if any(~isfinite(freq))
            incidencias(end+1,1) = sujeto + " | " + canal + " | Frecuencias no numericas";
            continue;
        end
        usar = (freq>=FMIN & freq<=FMAX);
        f = freq(usar);
        if isempty(f)
            incidencias(end+1,1) = sujeto + " | " + canal + " | Sin frecuencias 0.5-40 Hz";
            continue;
        end
        if numel(unique(f))~=numel(f) || any(diff(f)<=0)
            incidencias(end+1,1) = sujeto + " | " + canal + ...
                " | Frecuencias repetidas/no crecientes";
            continue;
        end

        nF = numel(f);
        X = nan(nF,N_MOM);
        for j=1:N_MOM
            ixCol = find(strcmpi(strtrim(nombres),NOM_MOMENTOS(j)),1);
            if isempty(ixCol)
                incidencias(end+1,1) = sujeto + " | " + canal + ...
                    " | Falta " + NOM_MOMENTOS(j);
            else
                v = convertirNumerico(T{:,ixCol});
                if numel(v)==numel(freq)
                    X(:,j)=v(usar);
                else
                    incidencias(end+1,1) = sujeto + " | " + canal + ...
                        " | Longitud incompatible de " + NOM_MOMENTOS(j);
                end
            end
        end

        % Una fila por sujeto x electrodo x momento x frecuencia.
        % Se conservan PSD NaN, ceros, negativas y valores extremos.
        nn = nF*N_MOM;
        momentosRep = repelem(CODIGOS_MOMENTOS,nF);
        fasesRep = repelem(FASES(:),nF);
        iniciosRep = repelem(INICIOS(:),nF);
        finesRep = repelem(FINES(:),nF);
        freqsRep = repmat(f,N_MOM,1);
        vals = X(:);
        valido = isfinite(vals) & vals>0;

        tablasSujeto{iE} = table( ...
            repmat(sujeto,nn,1),repmat(iE,nn,1), ...
            momentosRep,fasesRep,iniciosRep,finesRep,freqsRep,vals,valido, ...
            'VariableNames',{'Sujeto','Electrodo','Momento','Fase', ...
            'InicioOriginal_min','FinOriginal_min','Frecuencia_Hz', ...
            'PSD_lineal','PSD_FinitaPositiva'});

        sujetoQC(end+1:end+N_MOM,1) = repmat(sujeto,N_MOM,1);
        canalQC(end+1:end+N_MOM,1) = repmat(iE,N_MOM,1);
        momentoQC(end+1:end+N_MOM,1) = CODIGOS_MOMENTOS;
        nFrecuenciasQC(end+1:end+N_MOM,1) = nF;
        nValidosQC(end+1:end+N_MOM,1) = sum(isfinite(X) & X>0,1)';
        nInvalidosQC(end+1:end+N_MOM,1) = sum(isfinite(X) & X<=0,1)';
        nAusentesQC(end+1:end+N_MOM,1) = sum(~isfinite(X),1)';
    end

    tablasSujeto = tablasSujeto(~cellfun(@isempty,tablasSujeto));
    if isempty(tablasSujeto), continue; end
    Ts = vertcat(tablasSujeto{:});
    if primera
        writetable(Ts,OUT_LONG);
        primera = false;
    else
        writetable(Ts,OUT_LONG,'WriteMode','append','WriteVariableNames',false);
    end
    totalFilas = totalFilas + height(Ts);
    fprintf('  -> %d filas; acumuladas: %d\n',height(Ts),totalFilas);
end

if primera
    error('No se pudo importar ninguna PSD. Consultar archivos y nombres de hojas.');
end

QC = table(sujetoQC,canalQC,momentoQC,nFrecuenciasQC, ...
    nValidosQC,nInvalidosQC,nAusentesQC, ...
    'VariableNames',{'Sujeto','Electrodo','Momento','NumFrecuencias', ...
    'NumPSDPositivas','NumPSDNoPositivas','NumPSDNaNInf'});
writetable(QC,OUT_QC);

if ~isempty(incidencias)
    writetable(table(incidencias,'VariableNames',{'Incidencia'}),OUT_INCID);
    fprintf('ADVERTENCIA: %d incidencias. Ver: %s\n',numel(incidencias),OUT_INCID);
end

fprintf('\nLISTO: %d filas totales\n%s\n%s\n',totalFilas,OUT_LONG,OUT_QC);
fprintf('Equivalencias electrodos: %s\n',OUT_MAP);
fprintf('OC excluido. Codigos: OA=1; tSMS=2..16; POST=17..21.\n');
fprintf('No se han eliminado atipicos, agregado tiempos ni colapsado ROIs.\n');

%% FUNCION LOCAL
function v = convertirNumerico(x)
    if isnumeric(x) || islogical(x)
        v = double(x);
    else
        v = str2double(string(x));
    end
    v = v(:);
end
