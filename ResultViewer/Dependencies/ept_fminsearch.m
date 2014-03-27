function [x,fval,exitflag,output] = ept_fminsearch(funfcn,x,options,varargin)
%FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
%   X = FMINSEARCH(FUN,X0) starts at X0 and attempts to find a local minimizer 
%   X of the function FUN.  FUN is a function handle.  FUN accepts input X and 
%   returns a scalar function value F evaluated at X. X0 can be a scalar, vector 
%   or matrix.
%
%   X = FMINSEARCH(FUN,X0,OPTIONS)  minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created
%   with the OPTIMSET function.  See OPTIMSET for details.  FMINSEARCH uses
%   these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, FunValCheck,
%   PlotFcns, and OutputFcn.
%
%   X = FMINSEARCH(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fminsearch' in PROBLEM.solver. The PROBLEM structure must have
%   all the fields.
%
%   [X,FVAL]= FMINSEARCH(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINSEARCH(...) returns an EXITFLAG that describes 
%   the exit condition of FMINSEARCH. Possible values of EXITFLAG and the 
%   corresponding exit conditions are
%
%    1  Maximum coordinate difference between current best point and other
%       points in simplex is less than or equal to TolX, and corresponding 
%       difference in function values is less than or equal to TolFun.
%    0  Maximum number of function evaluations or iterations reached.
%   -1  Algorithm terminated by the output function.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINSEARCH(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name 
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fminsearch(@sin,3)
%     finds a minimum of the SIN function near 3.
%     In this case, SIN is a function that returns a scalar function value
%     SIN evaluated at X.
%
%     FUN can be an anonymous function:
%        X = fminsearch(@(x) norm(x),[1;2;3])
%     returns a point near the minimizer [0;0;0].
%
%     FUN can be a parameterized function. Use an anonymous function to
%     capture the problem-dependent parameters:
%        f = @(x,c) x(1).^2+c.*x(2).^2;  % The parameterized function.
%        c = 1.5;                        % The parameter.
%        X = fminsearch(@(x) f(x,c),[0.3;1])
%        
%   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.
%
%   See also OPTIMSET, FMINBND, FUNCTION_HANDLE.

%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
%   p.112-147, 1998.

%   Copyright 1984-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/02/14 03:33:29 $


defaultopt = struct('Display','notify','MaxIter','200*numberOfVariables',...
    'MaxFunEvals','200*numberOfVariables','TolX',1e-4,'TolFun',1e-4, ...
    'FunValCheck','off','OutputFcn',[],'PlotFcns',[]);

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(funfcn,'defaults')
    x = defaultopt;
    return
end

if nargin<3, options = []; end

% Detect problem structure input
if nargin == 1
    if isa(funfcn,'struct') 
        [funfcn,x,options] = separateOptimStruct(funfcn);
    else % Single input and non-structure
        error(message('MATLAB:fminsearch:InputArg'));
    end
end

if nargin == 0
    error(message('MATLAB:fminsearch:NotEnoughInputs'));
end


% Check for non-double inputs
if ~isa(x,'double')
  error(message('MATLAB:fminsearch:NonDoubleInput'))
end

n = numel(x);
numberOfVariables = n;

printtype = optimget(options,'Display',defaultopt,'fast');
tolx = optimget(options,'TolX',defaultopt,'fast');
tolf = optimget(options,'TolFun',defaultopt,'fast');
maxfun = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxiter = optimget(options,'MaxIter',defaultopt,'fast');
funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');

% In case the defaults were gathered from calling: optimset('fminsearch'):
if ischar(maxfun)
    if isequal(lower(maxfun),'200*numberofvariables')
        maxfun = 200*numberOfVariables;
    else
        error(message('MATLAB:fminsearch:OptMaxFunEvalsNotInteger'))
    end
end
if ischar(maxiter)
    if isequal(lower(maxiter),'200*numberofvariables')
        maxiter = 200*numberOfVariables;
    else
        error(message('MATLAB:fminsearch:OptMaxIterNotInteger'))
    end
end

switch printtype
    case {'notify','notify-detailed'}
        prnt = 1;
    case {'none','off'}
        prnt = 0;
    case {'iter','iter-detailed'}
        prnt = 3;
    case {'final','final-detailed'}
        prnt = 2;
    case 'simplex'
        prnt = 4;
    otherwise
        prnt = 1;
end
% Handle the output
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to plotfcns; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

header = ' Iteration   Func-count     min f(x)         Procedure';

% Convert to function handle as needed.
funfcn = fcnchk(funfcn,length(varargin));
% Add a wrapper function to check for Inf/NaN/complex values
if funValCheck
    % Add a wrapper function, CHECKFUN, to check for NaN/complex values without
    % having to change the calls that look like this:
    % f = funfcn(x,varargin{:});
    % x is the first argument to CHECKFUN, then the user's function,
    % then the elements of varargin. To accomplish this we need to add the 
    % user's function to the beginning of varargin, and change funfcn to be
    % CHECKFUN.
    varargin = {funfcn, varargin{:}};
    funfcn = @checkfun;
end

n = numel(x);

% Initialize parameters
rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
onesn = ones(1,n);
two2np1 = 2:n+1;
one2n = 1:n;

% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = zeros(n,n+1); fv = zeros(1,n+1);
v(:,1) = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
x(:) = xin;    % Change x to the form expected by funfcn
fv(:,1) = funfcn(x,varargin{:});
func_evals = 1;
itercount = 0;
how = '';
% Initial simplex setup continues later

% Initialize the output and plot functions.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'init',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end

% Print out initial f(x) as 0th iteration
if prnt == 3
    disp(' ')
    disp(header)
    fprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how);
elseif prnt == 4
    clc
    formatsave = get(0,{'format','formatspacing'});
    format compact
    format short e
    disp(' ')
    disp(how)
    v
    fv
    func_evals
end
% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop  % Stop per user request.
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end

% Continue setting up the initial simplex.
% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = 0.05;             % 5 percent deltas for non-zero terms
zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
for j = 1:n
    y = xin;
    if y(j) ~= 0
        y(j) = (1 + usual_delta)*y(j);
    else
        y(j) = zero_term_delta;
    end
    v(:,j+1) = y;
    x(:) = y; f = funfcn(x,varargin{:});
    fv(1,j+1) = f;
end

% sort so v(1,:) has the lowest function value
[fv,j] = sort(fv);
v = v(:,j);

how = 'initial simplex';
itercount = itercount + 1;
func_evals = n+1;
if prnt == 3
    fprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how)
elseif prnt == 4
    disp(' ')
    disp(how)
    v
    fv
    func_evals
end
% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop  % Stop per user request.
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end
exitflag = 1;

% Main algorithm: iterate until 
% (a) the maximum coordinate difference between the current best point and the 
% other points in the simplex is less than or equal to TolX. Specifically,
% until max(||v2-v1||,||v2-v1||,...,||v(n+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the 
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)
% The iteration stops if the maximum number of iterations or function evaluations 
% are exceeded
while func_evals < maxfun && itercount < maxiter
    if max(abs(fv(1)-fv(two2np1))) <= max(tolf,10*eps(fv(1))) && ...
            max(max(abs(v(:,two2np1)-v(:,onesn)))) <= max(tolx,10*eps(max(v(:,1))))
        break
    end
    
    % Compute the reflection point
    
    % xbar = average of the n (NOT n+1) best points
    xbar = sum(v(:,one2n), 2)/n;
    xr = (1 + rho)*xbar - rho*v(:,end);
    x(:) = xr; fxr = funfcn(x,varargin{:});
    func_evals = func_evals+1;
    
    if fxr < fv(:,1)
        % Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
        x(:) = xe; fxe = funfcn(x,varargin{:});
        func_evals = func_evals+1;
        if fxe < fxr
            v(:,end) = xe;
            fv(:,end) = fxe;
            how = 'expand';
        else
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        end
    else % fv(:,1) <= fxr
        if fxr < fv(:,n)
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        else % fxr >= fv(:,n)
            % Perform contraction
            if fxr < fv(:,end)
                % Perform an outside contraction
                xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
                x(:) = xc; fxc = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                
                if fxc <= fxr
                    v(:,end) = xc;
                    fv(:,end) = fxc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xcc = (1-psi)*xbar + psi*v(:,end);
                x(:) = xcc; fxcc = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                
                if fxcc < fv(:,end)
                    v(:,end) = xcc;
                    fv(:,end) = fxcc;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                for j=two2np1
                    v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
                    x(:) = v(:,j); fv(:,j) = funfcn(x,varargin{:});
                end
                func_evals = func_evals + n;
            end
        end
    end
    [fv,j] = sort(fv);
    v = v(:,j);
    itercount = itercount + 1;
    if prnt == 3
        fprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how)
    elseif prnt == 4
        disp(' ')
        disp(how)
        v
        fv
        func_evals
    end
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
            func_evals, how, fv(:,1),varargin{:});
        if stop  % Stop per user request.
            [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
            if  prnt > 0
                disp(output.message)
            end
            return;
        end
    end
end   % while

x(:) = v(:,1);
fval = fv(:,1);

if prnt == 4,
    % reset format
    set(0,{'format','formatspacing'},formatsave);
end
output.iterations = itercount;
output.funcCount = func_evals;
output.algorithm = 'Nelder-Mead simplex direct search';

% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'done',itercount, func_evals, how, fval, varargin{:});
end

if func_evals >= maxfun
    msg = getString(message('MATLAB:fminsearch:ExitingMaxFunctionEvals', sprintf('%f',fval)));
    if prnt > 0
        disp(' ')
        disp(msg)
    end
    exitflag = 0;
elseif itercount >= maxiter
    msg = getString(message('MATLAB:fminsearch:ExitingMaxIterations', sprintf('%f',fval)));
    if prnt > 0
        disp(' ')
        disp(msg)
    end
    exitflag = 0;
else
    msg = ...
      getString(message('MATLAB:fminsearch:OptimizationTerminatedXSatisfiesCriteria', ...
               sprintf('%e',tolx), sprintf('%e',tolf)));
    if prnt > 1
        disp(' ')
        disp(msg)
    end
    exitflag = 1;
end

output.message = msg;

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,state,iter,...
    numf,how,f,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.
%
% state - can have the values 'init','iter', or 'done'.

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
optimValues.iteration = iter;
optimValues.funccount = numf;
optimValues.fval = f;
optimValues.procedure = how;

xOutputfcn(:) = x;  % Set x to have user expected size
stop = false;
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error(message('MATLAB:fminsearch:InvalidState'))
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error(message('MATLAB:fminsearch:InvalidState'))
    end
end

%--------------------------------------------------------------------------
function [x,FVAL,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues)
% CLEANUPINTERRUPT updates or sets all the output arguments of FMINBND when the optimization
% is interrupted.

% Call plot function driver to finalize the plot function figure window. If
% no plot functions have been specified or the plot function figure no
% longer exists, this call just returns.
callAllOptimPlotFcns('cleanuponstopsignal');

x = xOutputfcn;
FVAL = optimValues.fval;
EXITFLAG = -1;
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.algorithm = 'Nelder-Mead simplex direct search';
OUTPUT.message = getString(message('MATLAB:fminsearch:OptimizationTerminatedPrematurelyByUser'));

%--------------------------------------------------------------------------
function f = checkfun(x,userfcn,varargin)
% CHECKFUN checks for complex or NaN results from userfcn.

f = userfcn(x,varargin{:});
% Note: we do not check for Inf as FMINSEARCH handles it naturally.
if isnan(f)
    error(message('MATLAB:fminsearch:checkfun:NaNFval', localChar( userfcn )));  
elseif ~isreal(f)
    error(message('MATLAB:fminsearch:checkfun:ComplexFval', localChar( userfcn )));  
end

%--------------------------------------------------------------------------
function strfcn = localChar(fcn)
% Convert the fcn to a string for printing

if ischar(fcn)
    strfcn = fcn;
elseif isa(fcn,'inline')
    strfcn = char(fcn);
elseif isa(fcn,'function_handle')
    strfcn = func2str(fcn);
else
    try
        strfcn = char(fcn);
    catch
        strfcn = getString(message('MATLAB:fminsearch:NameNotPrintable'));
    end
end


