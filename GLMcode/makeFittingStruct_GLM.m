function gg = makeFittingStruct_GLM(sta,DTsim,glmstruct,cellnumToFit);
% gg = makeFittingStruct_GLM(sta,DTsim,glmstruct,cellnumToFit);
%
% Initialize parameter structure for fitting of GLM model,
% normal parametrization of stim kernel
%
% Inputs:  sta = initial guess at kernel (or use all zeros if unknown)

% Set up structure
gg.k = [];
gg.dc = 0;
gg.ih = [];
gg.iht = [];
gg.ihbas = [];
gg.ihbas2 = [];
gg.ihbasprs = [];
gg.kt = [];
gg.ktbas = [];
gg.kbasprs = [];
gg.tsp = [];
gg.tspi = [];
gg.dt = DTsim;
gg.ih2 = [];
gg.ihbas2 = [];
gg.ihbasprs2 = [];
gg.tsp2 = [];
gg.couplednums = [];

% === Make temporal basis for stimulus filter =======================
[nkt,nkx] = size(sta);
% % ----- Set up temporal basis for stimulus kernel -----------
global KBASIS;
if isempty(KBASIS)
    KBASIS=10;
    fprintf('Number of K basis functions to use for kernel not declared, using default of 10\n')
end
global KFIRSTPEAK;
global KLASTPEAK;
if isempty(KFIRSTPEAK)
    KFIRSTPEAK=1;
    fprintf('Position of first K basis function peak undeclared, using default of frame 1\n');
end
if isempty(KLASTPEAK)
    KLASTPEAK=round(nkt*0.5);
    fprintf('Position of last K basis function peak undeclared, using default of (0.33*window length)\n');
end
global KLINEARITY;
if isempty(KLINEARITY)
    KLINEARITY=5;
    fprintf('KLinearity of peak distribution undeclared, using default of 5\n');
end
kbasprs.neye = 0; % Number of "identity" basis vectors near time of spike;
kbasprs.ncos = KBASIS; % Number of raised-cosine vectors to use  CHANGED FROM 3
kbasprs.kpeaks = [KFIRSTPEAK KLASTPEAK];  % Position of first and last bump (relative to identity bumps) CAHNGED FROM [10 ROUND(NKT*0.33)]
kbasprs.b =KLINEARITY; % Offset for nonlinear scaling (larger -> more linear)
ktbas = makeBasis_StimKernel(kbasprs,nkt);
gg.ktbas = ktbas;
gg.kbasprs = kbasprs;

% ======================================================================
% Set up basis for post-spike kernel
global HBASIS;
if isempty(HBASIS)
    HBASIS=10;
    fprintf('Number of I_h basis functions to use for kernel not declared, using default of 10\n')
end
global HFIRSTPEAK;
global HLASTPEAK;
if isempty(HFIRSTPEAK)
    HFIRSTPEAK=1;
    fprintf('Position of first I_h basis function peak undeclared, using default of frame 1\n');
end
if isempty(HLASTPEAK)
    HLASTPEAK=round(nkt*0.5);
    fprintf('Position of last I_h basis function peak undeclared, using default of (0.33*window length)\n');
end
global HLINEARITY;
if isempty(HLINEARITY)
    HLINEARITY=5;
    fprintf('HLinearity of peak distribution undeclared, using default of .5\n');
end
global ABSREF;
if isempty(ABSREF)
    ABSREF=1;
    fprintf('Absolute refractory period undeclared, using default of 1 frams\n');
end
ihbasprs.ncols = HBASIS;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [HFIRSTPEAK HLASTPEAK];  % Peak location for first and last vectors  [DTsim*10 2]
ihbasprs.b = HLINEARITY;  % How nonlinear to make spacings
ihbasprs.absref = ABSREF; %DTsim*10; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim);
gg.iht = iht;
gg.ihbas = ihbas;
gg.ihbasprs = ihbasprs;
gg.ih = zeros(size(ihbas,2),1);

% % ==================================================================
% set up initial K params
gg.kt = inv(gg.ktbas'*gg.ktbas)*gg.ktbas'*sta;
gg.k = gg.ktbas*gg.kt;

% % ==================================================================
% If full param struct passed in, match other params as well
if (nargin >= 3) 
    gg.dc = glmstruct.dc;
    gg.ih = glmstruct.ih;
    gg.iht = glmstruct.iht;

    %---Extract correct ih basis params, if present----
    if isfield(glmstruct, 'ihbasprs')
        if ~isempty(glmstruct.ihbasprs);
            ihbasprs = glmstruct.ihbasprs;
            [iht,ihbas] = makeBasis_PostSpike(ihbasprs,DTsim);

            % -- Do some error-checking ----
            if length(iht) ~= length(glmstruct.iht)
                error('mismatch between iht and h-kernel params ihbasprs');
            end
            if size(glmstruct.ih,2)>1 & (nargin < 4)
                error('multi-cell glm struct passed in without cell # to fit');
            end

            %--- Put glmstruct params into gg ----
            gg.iht = glmstruct.iht;
            gg.ihbas = ihbas;
            gg.ihbasprs = ihbasprs;
            if nargin == 3  % single-cell only
                gg.ih = inv(ihbas'*ihbas)*ihbas'*glmstruct.ih;
                gg.dc = glmstruct.dc;
            else % mulitcell-cell
                ncells = size(glmstruct.ih,2);
                ih0 = glmstruct.ih(:,:,cellnumToFit);
                ih1 = ih0(:,cellnumToFit);
                ih2 = ih0(:,setdiff(1:ncells,cellnumToFit));
                gg.ih = inv(ihbas'*ihbas)*ihbas'*[ih1 ih2];
                gg.dc = glmstruct.dc(cellnumToFit);
            end

        end
    end


end

