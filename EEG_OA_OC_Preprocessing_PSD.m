%% EEG OA/OC - PREPROCESAMIENTO Y PSD BASADO EN EEG_CREA v10
% Adaptacion al estudio actual:
% - 4 sujetos, sin grupos
% - Dos registros continuos por sujeto: OA y OC
% - Eliminar primeros y ultimos 5 segundos de cada EDF
% - 96 canales EEG; conservar F3
% - Sin ECG, fases, intervalos ni bloques
% - Filtrado -> deteccion de canales malos -> ASR -> ICA + ICLabel
%   -> interpolacion -> PSD Welch
% - ICLabel: umbral 0.70
% - Welch: ventana Hamming de 2 segundos y 50% de solapamiento
%
% IMPORTANTE:
% Ajustar EEGLAB_DIR y ROOT a las rutas reales.
% Usar el archivo de coordenadas correspondiente al montaje de 96 canales.
% Validar primero con un sujeto antes de procesar los cuatro.

clearvars;
clc;

%% 1. CONFIGURACION

EEGLAB_DIR = 'C:\TU_RUTA\eeglab';
ROOT = 'C:\TU_RUTA\DATOS_EEG';  % Contiene las carpetas de los 4 sujetos

DEST = fullfile(ROOT, 'EEG_preprocesados');
PSD_EXCEL_DIR = fullfile(DEST, 'PSD_por_sujeto');

% Dejar {} para procesar todas las carpetas de sujetos.
% Tambien puedes especificar: {'S01','S02','S03','S04'}
SUBJECT_FILTER = {};

EXPECTED_SUBJECTS = 4;
EXPECTED_EEG_CHANNELS = 96;

OVERWRITE = false;
SAVE_FULL_MAT = true;
SAVE_FILTERED_QC = true;
SAVE_PRE_ICA_QC = true;
SAVE_PSD_QC = true;
EXPORT_PSD_EXCEL = true;
VERIFY_SAVED_FILES = true;

% Recorte inicial y final
TRIM_SEC = 5;

% Tasa de muestreo esperada
FS_EXPECTED = 512;

% Filtro EEG
FILT_LOW = 0.5;
FILT_HIGH = 40;

% Parametros de limpieza del script original
FLATLINE = 5;
CHAN_CORR = 0.8;
LINE_NOISE = 4;
BURST_CRITERION = 20;
WINDOW_CRITERION = 0.25;
WINDOW_TOLERANCES = [-Inf 7];
ASR_DISTANCE = 'Euclidian';
ASR_MAXMEM = 64;
MIN_GOOD_CHANNELS = 8;

% ICLabel: conservar el umbral del script original
ICLABEL_THRESHOLD = 0.70;

% PSD mediante Welch
PSD_WINDOW_SEC = 2;
PSD_OVERLAP_FRACTION = 0.50;
PSD_QC_YLIM_COMMON = true;


%% 2. INICIALIZACION Y COMPROBACIONES

assert(isfile(fullfile(EEGLAB_DIR, 'eeglab.m')), ...
    'No se encontro eeglab.m en: %s', EEGLAB_DIR);

assert(isfolder(ROOT), ...
    'No existe la carpeta de datos: %s', ROOT);

addpath(EEGLAB_DIR);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui'); %#ok<ASGLU>

required = {'pop_biosig', 'pop_chanedit', 'pop_eegfiltnew', ...
    'pop_clean_rawdata', 'pop_runica', 'pop_iclabel', ...
    'pop_subcomp', 'pop_interp', 'pop_saveset', 'pop_select'};

for i = 1:numel(required)
    assert(~isempty(which(required{i})), ...
        'No esta disponible la funcion: %s', required{i});
end

% Archivo de localizaciones de referencia.
% Sustituir por las coordenadas reales del montaje si las tienes.
LOC = fullfile(fileparts(which('eeglab')), 'plugins', ...
    'dipfit', 'standard_BEM', 'elec', 'standard_1020.elc');

assert(isfile(LOC), ...
    'No se encuentra el archivo de localizaciones: %s', LOC);

if ~isfolder(DEST)
    mkdir(DEST);
end

if ~isfolder(PSD_EXCEL_DIR)
    mkdir(PSD_EXCEL_DIR);
end

runstamp = datestr(now, 'yyyymmdd_HHMMSS');
QCDIR = fullfile(DEST, ['AUDITORIA_' runstamp]);
mkdir(QCDIR);

diary(fullfile(QCDIR, 'CONSOLA_EJECUCION.txt'));
cleanupDiary = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('Inicio: %s\n', datestr(now));
fprintf('Origen: %s\n', ROOT);
fprintf('Destino: %s\n', DEST);
fprintf('EEGLAB: %s\n', which('eeglab'));

fprintf(['Filtro: %.1f-%.1f Hz | ICLabel: %.2f | ' ...
    'Welch: %.1f s, solapamiento %.0f%%\n'], ...
    FILT_LOW, FILT_HIGH, ICLABEL_THRESHOLD, ...
    PSD_WINDOW_SEC, 100*PSD_OVERLAP_FRACTION);

SUMMARY = table();
CHANNELS = table();
COMPONENTS = table();
INCIDENTS = table();
ALL_PSD_LONG = table();


%% 3. LOCALIZAR LOS SUJETOS (SIN GRUPOS)

subjectDirs = dir(ROOT);
subjectDirs = subjectDirs([subjectDirs.isdir] & ...
    ~ismember({subjectDirs.name}, {'.', '..'}));

% Excluir la carpeta de resultados
excludeDirs = {'EEG_preprocesados'};
subjectDirs = subjectDirs(~ismember({subjectDirs.name}, excludeDirs));

if ~isempty(SUBJECT_FILTER)
    subjectDirs = subjectDirs(ismember( ...
        lower(string({subjectDirs.name})), ...
        lower(string(SUBJECT_FILTER))));
end

if numel(subjectDirs) ~= EXPECTED_SUBJECTS
    warning('Se esperaban %d carpetas de sujetos y se encontraron %d.', ...
        EXPECTED_SUBJECTS, numel(subjectDirs));
end


%% 4. PROCESAR CADA SUJETO

for si = 1:numel(subjectDirs)

    subject = subjectDirs(si).name;
    srcSub = fullfile(ROOT, subject);

    EDFs = dir(fullfile(srcSub, '*.edf'));
    if isempty(EDFs)
        EDFs = dir(fullfile(srcSub, '*.EDF'));
    end

    if isempty(EDFs)
        warning('No se encontraron archivos EDF en %s', srcSub);
        INCIDENTS = [INCIDENTS; table(string(subject), "", ...
            "No hay archivos EDF", ...
            'VariableNames', {'Sujeto','Condicion','Incidencia'})]; %#ok<AGROW>
        continue;
    end

    outSub = fullfile(DEST, subject);
    workSub = fullfile(outSub, 'INTERMEDIOS');

    if ~isfolder(workSub)
        mkdir(workSub);
    end

    % Identificar un EDF OA y un EDF OC por el nombre del archivo
    fileOA = '';
    fileOC = '';

    for fi = 1:numel(EDFs)

        [~, stemTmp] = fileparts(EDFs(fi).name);
        nameU = upper(stemTmp);

        isOA = ~isempty(regexp(nameU, ...
            '(^|[^A-Z0-9])OA([^A-Z0-9]|$)', 'once'));

        isOC = ~isempty(regexp(nameU, ...
            '(^|[^A-Z0-9])OC([^A-Z0-9]|$)', 'once'));

        if isOA && isOC
            error('El nombre contiene OA y OC a la vez: %s', ...
                EDFs(fi).name);

        elseif isOA
            assert(isempty(fileOA), ...
                'Hay mas de un EDF OA en %s', srcSub);
            fileOA = fullfile(srcSub, EDFs(fi).name);

        elseif isOC
            assert(isempty(fileOC), ...
                'Hay mas de un EDF OC en %s', srcSub);
            fileOC = fullfile(srcSub, EDFs(fi).name);
        end
    end

    if isempty(fileOA) || isempty(fileOC)
        warning('Falta OA o OC para %s. Se omite este sujeto.', subject);

        INCIDENTS = [INCIDENTS; table(string(subject), "", ...
            "Falta OA o OC", ...
            'VariableNames', {'Sujeto','Condicion','Incidencia'})]; %#ok<AGROW>
        continue;
    end

    fprintf('\n============================================\n');
    fprintf('SUJETO: %s\n', subject);
    fprintf('============================================\n');

    conditions = {'OA', 'OC'};
    srcFiles = {fileOA, fileOC};
    PSD_COND = struct();

    %% 4.1 PROCESAR OA Y OC

    for ci = 1:2

        condition = conditions{ci};
        srcEDF = srcFiles{ci};
        [~, stem] = fileparts(srcEDF);

        outfile = fullfile(workSub, [stem '_preprocesado.set']);
        psdFile = fullfile(workSub, [stem '_PSD_EXPORT.mat']);

        fprintf('\n===== %s | %s =====\n', subject, condition);

        if isfile(outfile) && isfile(psdFile) && ~OVERWRITE
            fprintf('Resultados existentes; se cargan: %s\n', stem);

            tmp = load(psdFile, 'P');
            PSD_COND.(condition) = tmp.P;

            Tlong = psd_to_long_table(tmp.P, subject, condition);
            ALL_PSD_LONG = [ALL_PSD_LONG; Tlong]; %#ok<AGROW>
            continue;
        end

        try

            %% 4.1.1 IMPORTAR EDF

            EEG = pop_biosig(srcEDF);
            EEG = eeg_checkset(EEG);

            assert(EEG.trials == 1, ...
                'Se esperaba un registro continuo; trials=%d', EEG.trials);

            assert(abs(EEG.srate-FS_EXPECTED) < 1e-6, ...
                ['Tasa de muestreo inesperada: %.3f Hz. ' ...
                 'Se esperaba %.1f Hz.'], EEG.srate, FS_EXPECTED);

            fs = EEG.srate;

            %% 4.1.2 EXCLUIR CANALES AUXILIARES

            rawLabels = string({EEG.chanlocs.labels});

            % No se realiza procesamiento ECG.
            % Se excluyen solo canales auxiliares identificables
            % por sus etiquetas.
            auxIdx = find( ...
                contains(upper(rawLabels), 'MKR') | ...
                contains(upper(rawLabels), 'MARKER') | ...
                contains(upper(rawLabels), 'TRIGGER') | ...
                strcmpi(rawLabels, 'STATUS') | ...
                strcmpi(rawLabels, 'TRIG'));

            if ~isempty(auxIdx)
                fprintf('Canales auxiliares excluidos: %s\n', ...
                    strjoin(rawLabels(auxIdx), ', '));

                EEG = pop_select(EEG, ...
                    'nochannel', cellstr(rawLabels(auxIdx)));
                EEG = eeg_checkset(EEG);
            end

            assert(EEG.nbchan == EXPECTED_EEG_CHANNELS, ...
                ['Se esperaban %d canales EEG tras excluir auxiliares, ' ...
                 'pero se encontraron %d. Revisa las etiquetas.'], ...
                EXPECTED_EEG_CHANNELS, EEG.nbchan);


            %% 4.1.3 ELIMINAR LOS PRIMEROS Y ULTIMOS 5 SEGUNDOS

            Noriginal = EEG.pnts;

            assert(Noriginal > 2*round(TRIM_SEC*fs), ...
                'Registro demasiado corto para aplicar el recorte.');

            firstSample = round(TRIM_SEC*fs) + 1;
            lastSample = Noriginal - round(TRIM_SEC*fs);

            EEG = pop_select(EEG, ...
                'point', [firstSample lastSample]);
            EEG = eeg_checkset(EEG);

            N = EEG.pnts;
            trimmedDuration = N/fs;

            fprintf('Duracion original: %.2f s\n', Noriginal/fs);
            fprintf('Duracion tras recorte: %.2f s\n', trimmedDuration);

            assert(abs(trimmedDuration - ...
                (Noriginal/fs - 2*TRIM_SEC)) <= 2/fs, ...
                'Duracion inesperada despues del recorte.');


            %% 4.1.4 LOCALIZACIONES DE ELECTRODOS

            EEG = pop_chanedit(EEG, 'lookup', LOC);
            EEG = eeg_checkset(EEG);

            referenceChanlocs = EEG.chanlocs;
            originalLabels = string({EEG.chanlocs.labels});

            hasXYZ = arrayfun(@(x) ...
                isfield(x,'X') && isfield(x,'Y') && isfield(x,'Z') && ...
                ~isempty(x.X) && ~isempty(x.Y) && ~isempty(x.Z), ...
                EEG.chanlocs);

            if ~all(hasXYZ)
                missingLocs = originalLabels(~hasXYZ);

                warning(['Algunos canales no tienen coordenadas: %s. ' ...
                    'Comprueba el archivo de localizaciones del montaje.'], ...
                    strjoin(missingLocs, ', '));
            end


            %% 4.1.5 FILTRADO 0.5-40 Hz

            EEG = pop_eegfiltnew(EEG, ...
                'locutoff', FILT_LOW, 'plotfreqz', 0);

            EEG = pop_eegfiltnew(EEG, ...
                'hicutoff', FILT_HIGH, 'plotfreqz', 0);

            assert(EEG.pnts == N, ...
                'El filtrado ha cambiado la longitud de la señal.');

            filtered = EEG;

            if SAVE_FILTERED_QC
                qc = filtered;
                qc.setname = [stem '_filtrado'];

                pop_saveset(qc, ...
                    'filename', [stem '_filtrado.set'], ...
                    'filepath', workSub);
            end


            %% 4.1.6 DETECTAR CANALES MALOS SIN ASR

            detect = pop_clean_rawdata(filtered, ...
                'FlatlineCriterion', FLATLINE, ...
                'ChannelCriterion', CHAN_CORR, ...
                'LineNoiseCriterion', LINE_NOISE, ...
                'Highpass', 'off', ...
                'BurstCriterion', 'off', ...
                'WindowCriterion', 'off');

            goodLabels = string({detect.chanlocs.labels});
            badLabels = setdiff(originalLabels, goodLabels, 'stable');

            assert(numel(goodLabels) >= MIN_GOOD_CHANNELS, ...
                'Solo %d canales validos; revisa la limpieza.', ...
                numel(goodLabels));

            assert(detect.pnts == N, ...
                'La deteccion de canales ha alterado las muestras.');

            % Excluir los canales malos de los datos filtrados originales
            EEG = pop_select(filtered, ...
                'channel', cellstr(goodLabels));

            assert(EEG.pnts == N, ...
                'La exclusion de canales ha cambiado las muestras.');

            for k = 1:numel(originalLabels)
                CHANNELS = [CHANNELS; table( ...
                    string(subject), string(condition), originalLabels(k), ...
                    ismember(originalLabels(k), badLabels), ...
                    'VariableNames', ...
                    {'Sujeto','Condicion','Canal','EliminadoAntesASR'})]; %#ok<AGROW>
            end


            %% 4.1.7 ASR CON RECHAZO

            EEG = pop_clean_rawdata(EEG, ...
                'FlatlineCriterion', 'off', ...
                'ChannelCriterion', 'off', ...
                'LineNoiseCriterion', 'off', ...
                'Highpass', 'off', ...
                'BurstCriterion', BURST_CRITERION, ...
                'WindowCriterion', WINDOW_CRITERION, ...
                'BurstRejection', 'on', ...
                'Distance', ASR_DISTANCE, ...
                'WindowCriterionTolerances', WINDOW_TOLERANCES, ...
                'MaxMem', ASR_MAXMEM);

            assert(isfield(EEG.etc, 'clean_sample_mask'), ...
                'ASR no ha generado clean_sample_mask.');

            mask = logical(EEG.etc.clean_sample_mask(:)');

            assert(numel(mask) == N, ...
                'La mascara ASR no coincide con las muestras originales.');

            assert(nnz(mask) == EEG.pnts, ...
                'La mascara ASR no coincide con EEG.pnts.');

            kept_orig_idx = find(mask);

            assert(all(diff(kept_orig_idx) > 0), ...
                'Los indices conservados por ASR no son crecientes.');

            if SAVE_PRE_ICA_QC
                qc = EEG;
                qc.setname = [stem '_postASR_preICA'];

                pop_saveset(qc, ...
                    'filename', [stem '_postASR_preICA.set'], ...
                    'filepath', workSub);
            end


            %% 4.1.8 ICA

            Xrank = double(EEG.data(:, ...
                1:min(EEG.pnts, round(30*fs))));

            datarank = rank(Xrank*Xrank');

            assert(datarank >= 2, ...
                'Rango insuficiente para realizar ICA.');

            if datarank < EEG.nbchan
                fprintf('Rango ICA %d < %d canales; se usa PCA=%d.\n', ...
                    datarank, EEG.nbchan, datarank);

                EEG = pop_runica(EEG, ...
                    'icatype', 'runica', ...
                    'extended', 1, ...
                    'pca', datarank, ...
                    'interrupt', 'off');
            else
                EEG = pop_runica(EEG, ...
                    'icatype', 'runica', ...
                    'extended', 1, ...
                    'interrupt', 'off');
            end


            %% 4.1.9 CLASIFICACION ICLabel

            EEG = pop_iclabel(EEG, 'default');

            probs = EEG.etc.ic_classification.ICLabel.classifications;
            classes = string(EEG.etc.ic_classification.ICLabel.classes);

            assert(size(probs,2) >= 7, ...
                'ICLabel no ha devuelto las siete clases.');

            assert(size(probs,1) == size(EEG.icaweights,1), ...
                'El numero de componentes IC no coincide.');


            %% 4.1.10 ELIMINAR COMPONENTES ARTEFACTUALES

            % Clases ICLabel:
            % 1 Brain, 2 Muscle, 3 Eye, 4 Heart,
            % 5 Line Noise, 6 Channel Noise, 7 Other.
            %
            % Se conserva el umbral 0.70 del script original.
            % No se elimina la clase "Other" de forma automatica.

            artifactClasses = 2:6;
            toRemove = any( ...
                probs(:, artifactClasses) >= ICLABEL_THRESHOLD, 2);

            removedICs = find(toRemove);

            for ic = 1:size(probs,1)

                [mx, imx] = max(probs(ic,:));

                COMPONENTS = [COMPONENTS; table( ...
                    string(subject), string(condition), ic, ...
                    classes(imx), mx, ...
                    probs(ic,1), probs(ic,2), probs(ic,3), ...
                    probs(ic,4), probs(ic,5), probs(ic,6), probs(ic,7), ...
                    toRemove(ic), ...
                    'VariableNames', ...
                    {'Sujeto','Condicion','IC','ClasePrincipal','ProbPrincipal', ...
                    'PBrain','PMuscle','PEye','PHeart','PLineNoise', ...
                    'PChannelNoise','POther','Eliminada'})]; %#ok<AGROW>
            end

            assert(nnz(~toRemove) >= 2, ...
                ['ICLabel marcaria casi todas las IC como artefactos. ' ...
                 'Revisar manualmente antes de continuar.']);

            if ~isempty(removedICs)
                EEG = pop_subcomp(EEG, removedICs, 0);
            end


            %% 4.1.11 INTERPOLAR LOS CANALES MALOS

            if ~isempty(badLabels)
                EEG = pop_interp(EEG, referenceChanlocs, 'spherical');
            end

            EEG = eeg_checkset(EEG);

            assert(EEG.nbchan == numel(referenceChanlocs), ...
                'No se ha restaurado el montaje completo tras interpolar.');

            assert(EEG.pnts == nnz(mask), ...
                'ICA/interpolacion han alterado las muestras conservadas.');

            assert(any(strcmpi(string({EEG.chanlocs.labels}), 'F3')), ...
                'F3 no aparece en el resultado final.');


            %% 4.1.12 GUARDAR EL EEG PREPROCESADO

            EEG.setname = [stem '_preprocesado'];
            EEG.etc.clean_sample_mask = mask;

            EEG.etc.EEG_OAOC = struct( ...
                'version', 'OA_OC_96ch_v1', ...
                'sourceEDF', srcEDF, ...
                'subject', subject, ...
                'condition', condition, ...
                'originalSamplesAfterTrim', N, ...
                'originalSamplesBeforeTrim', Noriginal, ...
                'sampleRate', fs, ...
                'trimSecondsAtEachEnd', TRIM_SEC, ...
                'analyzedDuration_s', trimmedDuration, ...
                'badChannels', {cellstr(badLabels)}, ...
                'originalChanlocs', {referenceChanlocs}, ...
                'originalSampleIndexAfterASR', kept_orig_idx, ...
                'parameters', struct( ...
                    'filterHz', [FILT_LOW FILT_HIGH], ...
                    'FlatlineCriterion', FLATLINE, ...
                    'ChannelCriterion', CHAN_CORR, ...
                    'LineNoiseCriterion', LINE_NOISE, ...
                    'BurstCriterion', BURST_CRITERION, ...
                    'WindowCriterion', WINDOW_CRITERION, ...
                    'WindowCriterionTolerances', WINDOW_TOLERANCES, ...
                    'BurstRejection', 'on', ...
                    'ICLabelThreshold', ICLABEL_THRESHOLD, ...
                    'PSDwindow_s', PSD_WINDOW_SEC, ...
                    'PSDoverlap', PSD_OVERLAP_FRACTION));

            maskOrig = mask; %#ok<NASGU>
            original_sample_index = kept_orig_idx; %#ok<NASGU>
            original_srate = fs; %#ok<NASGU>
            original_n_samples = N; %#ok<NASGU>

            save(fullfile(workSub, [stem '_MAPA_TEMPORAL.mat']), ...
                'maskOrig', 'original_sample_index', ...
                'original_srate', 'original_n_samples', '-v7.3');

            pop_saveset(EEG, ...
                'filename', [stem '_preprocesado.set'], ...
                'filepath', workSub);

            if SAVE_FULL_MAT
                save(fullfile(outSub, ['EEG_preprocess_' stem '.mat']), ...
                    'EEG', '-v7.3');
            end

            if VERIFY_SAVED_FILES
                R = pop_loadset( ...
                    'filename', [stem '_preprocesado.set'], ...
                    'filepath', workSub);

                assert(R.nbchan == EEG.nbchan && R.pnts == EEG.pnts, ...
                    'El .set recargado no coincide con el original.');

                assert(isequal(string({R.chanlocs.labels}), ...
                    string({EEG.chanlocs.labels})), ...
                    'Las etiquetas cambiaron al recargar el .set.');
            end


            %% 4.1.13 PSD DE CONTROL ANTES/DESPUES

            if SAVE_PSD_QC
                try
                    qcinfo = plot_psd_qc_oaoc( ...
                        filtered, EEG, mask, referenceChanlocs, ...
                        fs, stem, outSub, badLabels, ...
                        PSD_WINDOW_SEC, PSD_OVERLAP_FRACTION, ...
                        PSD_QC_YLIM_COMMON);

                    fprintf('PSD QC guardada: %s\n', qcinfo.file);

                catch ME_QC
                    warning('Fallo la figura PSD QC de %s: %s', ...
                        stem, ME_QC.message);

                    INCIDENTS = [INCIDENTS; table( ...
                        string(subject), string(condition), ...
                        string(['PSD_QC: ' ME_QC.message]), ...
                        'VariableNames', ...
                        {'Sujeto','Condicion','Incidencia'})]; %#ok<AGROW>
                end
            end


            %% 4.1.14 CALCULAR PSD PARA OA U OC

            P = calculate_psd_oaoc( ...
                EEG, mask, stem, originalLabels, fs, ...
                PSD_WINDOW_SEC, PSD_OVERLAP_FRACTION);

            P.subject = subject;
            P.condition = condition;
            P.sourceEDF = srcEDF;
            P.trimSecondsAtEachEnd = TRIM_SEC;
            P.originalDuration_s = Noriginal/fs;
            P.analyzedDuration_s = trimmedDuration;
            P.badChannels = cellstr(badLabels);

            save(psdFile, 'P', '-v7.3');

            PSD_COND.(condition) = P;

            Tlong = psd_to_long_table(P, subject, condition);
            ALL_PSD_LONG = [ALL_PSD_LONG; Tlong]; %#ok<AGROW>

            SUMMARY = [SUMMARY; table( ...
                string(subject), string(condition), ...
                Noriginal/fs, trimmedDuration, ...
                nnz(mask), N-nnz(mask), ...
                100*(1-nnz(mask)/N), ...
                numel(originalLabels), numel(goodLabels), ...
                numel(badLabels), string(strjoin(badLabels, ',')), ...
                size(probs,1), nnz(toRemove), ...
                string(mat2str(removedICs(:)')), P.nWelchWindows, ...
                'VariableNames', ...
                {'Sujeto','Condicion','DuracionOriginal_s', ...
                'DuracionAnalizada_s','MuestrasConservadasASR', ...
                'MuestrasRechazadasASR','RechazoASR_pct', ...
                'CanalesOriginales','CanalesParaASR_ICA', ...
                'NumCanalesMalos','CanalesMalos','NumIC', ...
                'NumICEliminadas','ICEliminadas','VentanasWelch'})]; %#ok<AGROW>

            fprintf(['OK %s %s: rechazo ASR %.2f%%; canales malos {%s}; ' ...
                'IC eliminadas %d/%d; ventanas Welch %d.\n'], ...
                subject, condition, ...
                100*(1-nnz(mask)/N), strjoin(badLabels, ','), ...
                nnz(toRemove), size(probs,1), P.nWelchWindows);

        catch ME
            fprintf(2, 'ERROR en %s (%s): %s\n', ...
                subject, condition, ME.getReport('basic'));

            INCIDENTS = [INCIDENTS; table( ...
                string(subject), string(condition), string(ME.message), ...
                'VariableNames', ...
                {'Sujeto','Condicion','Incidencia'})]; %#ok<AGROW>
        end
    end


    %% 4.2 EXPORTAR EXCEL OA VS OC PARA CADA SUJETO

    if EXPORT_PSD_EXCEL && ...
            isfield(PSD_COND, 'OA') && isfield(PSD_COND, 'OC')

        P_OA = PSD_COND.OA;
        P_OC = PSD_COND.OC;

        export_subject_oaoc( ...
            P_OA, P_OC, PSD_EXCEL_DIR, subject);

        fprintf('Excel OA/OC exportado para %s.\n', subject);
    end
end


%% 5. EXPORTAR TABLAS GLOBALES Y AUDITORIA

if ~isempty(ALL_PSD_LONG)
    writetable(ALL_PSD_LONG, ...
        fullfile(PSD_EXCEL_DIR, 'PSD_TODOS_OA_OC_long.csv'));
end

if ~isempty(SUMMARY)
    writetable(SUMMARY, fullfile(QCDIR, 'RESUMEN.csv'));
end

if ~isempty(CHANNELS)
    writetable(CHANNELS, fullfile(QCDIR, 'CANALES.csv'));
end

if ~isempty(COMPONENTS)
    writetable(COMPONENTS, ...
        fullfile(QCDIR, 'COMPONENTES_ICLABEL.csv'));
end

if ~isempty(INCIDENTS)
    writetable(INCIDENTS, fullfile(QCDIR, 'INCIDENCIAS.csv'));
end

fprintf('\nTERMINADO.\nResultados: %s\nAuditoria: %s\n', DEST, QCDIR);


%% ================================================================
% FUNCIONES LOCALES
%% ================================================================

function P = calculate_psd_oaoc(EEG, mask, stem, ...
    originalLabels, fs, welchSec, overlapFraction)

    % PSD Welch calculada dentro de cada tramo continuo conservado por ASR.
    % No se crean intervalos experimentales y ninguna ventana cruza huecos.

    labels = string({EEG.chanlocs.labels});
    [present, loc] = ismember(string(originalLabels), labels);

    assert(all(present), ...
        'PSD: faltan canales tras la interpolacion.');

    X = double(EEG.data(loc,:));
    kept = find(logical(mask(:)'));

    assert(numel(kept) == size(X,2), ...
        'PSD: mascara y muestras no coinciden.');

    L = round(fs*welchSec);
    hop = round(L*(1-overlapFraction));
    nfft = L;
    win = hamming(L);

    freq = (0:nfft/2)'*fs/nfft;
    nCh = numel(originalLabels);

    accum = zeros(nCh, numel(freq));
    nTotalWindows = 0;
    fragmentCount = 0;

    if isempty(kept)
        error('PSD: ASR no ha conservado muestras.');
    end

    starts = [1 find(diff(kept)>1)+1];
    ends = [starts(2:end)-1 numel(kept)];

    for r = 1:numel(starts)

        idx = starts(r):ends(r);

        if numel(idx) < L
            continue;
        end

        fragmentCount = fragmentCount + 1;
        nWin = floor((numel(idx)-L)/hop)+1;

        for ch = 1:nCh

            [pxx, fActual] = pwelch( ...
                X(ch,idx), win, L-hop, nfft, fs);

            assert(numel(fActual) == numel(freq) && ...
                max(abs(fActual(:)-freq)) < 1e-10, ...
                'PSD: eje de frecuencias inesperado.');

            accum(ch,:) = accum(ch,:) + nWin*pxx(:)';
        end

        nTotalWindows = nTotalWindows+nWin;
    end

    if nTotalWindows > 0
        psd = accum/nTotalWindows;
    else
        psd = nan(nCh, numel(freq));
    end

    P = struct();
    P.stem = stem;
    P.channelLabels = cellstr(string(originalLabels));
    P.frequencies = freq;
    P.psd = psd;
    P.fs = fs;
    P.nSamplesOriginalTrimmed = numel(mask);
    P.nSamplesKeptASR = nnz(mask);
    P.nSamplesRejectedASR = numel(mask)-nnz(mask);
    P.coverageASR_pct = 100*nnz(mask)/numel(mask);
    P.nWelchWindows = nTotalWindows;
    P.nContinuousFragmentsWithPSD = fragmentCount;
    P.welchWindowSeconds = welchSec;
    P.overlapFraction = overlapFraction;
end


function T = psd_to_long_table(P, subject, condition)

    nCh = size(P.psd,1);
    nF = numel(P.frequencies);

    channels = repelem(string(P.channelLabels(:)), nF);
    frequencies = repmat(P.frequencies(:), nCh, 1);
    psdVals = reshape(P.psd', [], 1);

    T = table( ...
        repmat(string(subject), nCh*nF, 1), ...
        repmat(string(condition), nCh*nF, 1), ...
        channels, frequencies, psdVals, ...
        'VariableNames', ...
        {'Sujeto','Condicion','Canal','Frecuencia_Hz','PSD'});
end


function export_subject_oaoc(P_OA, P_OC, folder, subject)

    assert(isequal(string(P_OA.channelLabels), ...
        string(P_OC.channelLabels)), ...
        'El montaje de OA y OC no coincide.');

    assert(isequal(P_OA.frequencies, P_OC.frequencies), ...
        'Las frecuencias de OA y OC no coinciden.');

    if ~isfolder(folder)
        mkdir(folder);
    end

    xlsx = fullfile(folder, [subject '_PSD_OA_OC.xlsx']);

    if isfile(xlsx)
        delete(xlsx);
    end

    % Una hoja por canal:
    % Frecuencia_Hz, PSD_OA y PSD_OC.
    for ch = 1:numel(P_OA.channelLabels)

        label = matlab.lang.makeValidName( ...
            char(P_OA.channelLabels{ch}));

        if strlength(label) > 31
            label = extractBetween(string(label), 1, 31);
            label = char(label);
        end

        T = table( ...
            P_OA.frequencies(:), ...
            P_OA.psd(ch,:)', ...
            P_OC.psd(ch,:)', ...
            'VariableNames', ...
            {'Frecuencia_Hz','PSD_OA','PSD_OC'});

        writetable(T, xlsx, 'Sheet', label);
    end
end


function info = plot_psd_qc_oaoc(filtered, cleaned, mask, ...
    allChanlocs, fs, stem, dest, badLabels, ...
    welchSec, overlapFraction, commonY)

    labels = string({allChanlocs.labels});

    [okBefore, ordBefore] = ismember( ...
        labels, string({filtered.chanlocs.labels}));

    [okAfter, ordAfter] = ismember( ...
        labels, string({cleaned.chanlocs.labels}));

    assert(all(okBefore) && all(okAfter), ...
        'PSD QC: faltan canales.');

    Xbefore = double(filtered.data(ordBefore,:));
    Xafter = double(cleaned.data(ordAfter,:));

    L = round(welchSec*fs);
    hop = round(L*(1-overlapFraction));
    win = hamming(L);
    nfft = L;

    [~, f] = pwelch(Xbefore(1,:), win, L-hop, nfft, fs);

    psdBefore = nan(numel(labels), numel(f));

    for ch = 1:numel(labels)
        [pBeforeCh, ~] = pwelch( ...
            Xbefore(ch,:), win, L-hop, nfft, fs);
        psdBefore(ch,:) = pBeforeCh(:)';
    end

    kept = find(logical(mask(:)'));
    starts = [1 find(diff(kept)>1)+1];
    ends = [starts(2:end)-1 numel(kept)];

    accum = zeros(numel(labels), numel(f));
    nWinTotal = 0;

    for r = 1:numel(starts)

        idx = starts(r):ends(r);

        if numel(idx) < L
            continue;
        end

        nWin = floor((numel(idx)-L)/hop)+1;

        for ch = 1:numel(labels)

            [p, ff] = pwelch( ...
                Xafter(ch,idx), win, L-hop, nfft, fs);

            assert(max(abs(ff(:)-f(:))) < 1e-10, ...
                'PSD QC: frecuencias incompatibles.');

            accum(ch,:) = accum(ch,:) + nWin*p(:)';
        end

        nWinTotal = nWinTotal+nWin;
    end

    if nWinTotal > 0
        psdAfter = accum/nWinTotal;
    else
        psdAfter = nan(size(accum));
    end

    use = f >= 0.5 & f <= 40;

    dbBefore = 10*log10(max(psdBefore, realmin));
    dbAfter = 10*log10(max(psdAfter, realmin));

    colors = lines(numel(labels));

    h = figure('Visible','off','Color','w', ...
        'Position',[60 80 1600 690]);

    cleanupFig = onCleanup(@() close(h)); %#ok<NASGU>

    ax1 = subplot(1,2,1,'Parent',h);
    hold(ax1,'on');

    ax2 = subplot(1,2,2,'Parent',h);
    hold(ax2,'on');

    for ch = 1:numel(labels)
        plot(ax1, f(use), dbBefore(ch,use), ...
            'LineWidth',0.8,'Color',colors(ch,:));

        plot(ax2, f(use), dbAfter(ch,use), ...
            'LineWidth',0.8,'Color',colors(ch,:));
    end

    title(ax1, 'Antes: filtrado 0.5-40 Hz');
    title(ax2, sprintf( ...
        'Despues: ASR + ICA + interpolacion (%d ventanas Welch)', ...
        nWinTotal));

    for ax = [ax1 ax2]
        xlabel(ax, 'Frecuencia (Hz)');
        ylabel(ax, 'PSD (dB/Hz)');
        xlim(ax, [0.5 40]);
        grid(ax, 'on');
        box(ax, 'on');
    end

    if commonY
        vals = [dbBefore(:,use), dbAfter(:,use)];
        vals = vals(isfinite(vals));

        if ~isempty(vals)
            yl = [floor(min(vals)) ceil(max(vals))];

            if yl(1) == yl(2)
                yl = yl + [-1 1];
            end

            ylim(ax1, yl);
            ylim(ax2, yl);
        end
    end

    legend(ax2, cellstr(labels), ...
        'Location','eastoutside','FontSize',7);

    annotation(h, 'textbox', [0.02 0.94 0.96 0.05], ...
        'String', sprintf( ...
        '%s | PSD QC antes/despues | canales malos: %s', ...
        stem, strjoin(badLabels, ', ')), ...
        'LineStyle','none','FontSize',12,'FontWeight','bold', ...
        'HorizontalAlignment','center','Interpreter','none');

    annotation(h, 'textbox', [0.02 0.005 0.96 0.055], ...
        'String', sprintf( ...
        ['Welch: Hamming %.1f s, solapamiento %.0f%%. ' ...
         'PSD posterior calculada dentro de tramos continuos conservados por ASR.'], ...
        welchSec, 100*overlapFraction), ...
        'LineStyle','none','FontSize',9, ...
        'HorizontalAlignment','center');

    file = fullfile(dest, ...
        [stem '_PSD_QC_ANTES_DESPUES.png']);

    print(h, file, '-dpng', '-r150');

    info = struct('file',file,'nWindowsAfter',nWinTotal);
end