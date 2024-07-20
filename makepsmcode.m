function [newfunc, valuelist, fvalue, location] ...
    = makepsmcode(fhandle, numDEs, t0, y0)
% Evaluates fhandle and writes a .m file to generate series coefficients
%    numDEs: the number of functions (number of scalar ODEs) in fhandle
%    t0, y0: OPTIONAL initial values if needed to avoid domain error;
%            if omitted, t0,y0 are randomly generated in interval (0,1)
% Output: program to generate series coefficients for any point and order
%    .m file with newfuncname being fhandle name augmented by 'series'
% Other (optional) output arguments
%    newfuncname: string for function name of file created or overwritten
%    valuelist: all intermediate values computed in f(t0,y0) evaluation.
%    value: f(t0,y0) final value, may be a vector.
%    location: index of final value(s) within the valuelist.
% by Richard Neidinger 1/8/21; rev 10/26/22(allowed DE = constant);
% 6/8/24 for better code presentation, finalized for publication 7/19/24.

% if t0 and y0 are not given, we can use any random values
if nargin == 2  % need some t0,y0 value to go through evaluation steps
    t0 = rand;  % random in (0,1)
    y0 = rand(numDEs,1);
end
clear trace  % restart trace class, clears all persistent variables
t = trace(t0, 0);  % 0 indicates independent variable
y = trace(y0, 1);  % 1 indicates dependent variable in DE for y(t)

h = feval(fhandle, t, y);  % In evaluation of trace objects,  
      % each overloaded operation writes string(s) of MATLAB code for the
      % operation into oplist and corresponding recurrence in recurlist.
if ~isa(h,'trace')
    disp(['Overloading ',func2str(fhandle),' FAILED, try omitting preallocation.']);
end
[fvalue, location] = valueandopcount(h);  % vectors for f(t,y) results
[numops, valuelist, oplist, recurlist] = trace.count_val;%#ok<ASGLU>
% %#ok<ASGLU> surpresses a Code Analyzer warning about unused assignments. 

% Augment cell arrays oplist and recurlist with code to update Y by the DE.
oplist = [oplist; {'% Update Y linear coefficient by y'' = f(t,y)'}];
recurlist = [recurlist; {'% Update Y next coefficient using y'' = f(t,y)'}];
for i = 1:numDEs
    if isnan(location(i)) % trivial case of DE = constant
        derivname = num2str(fvalue(i));
        seriesname = '0';
    else
        derivname = [trace.name(location(i)),',1)'];
        seriesname = [trace.name(location(i)),',j)/j'];
    end
    yprime = ['Y(',num2str(i),',2) = ',derivname,';'];
    oplist = [oplist; {yprime}];
    nextterm = ['Y(',num2str(i),',j+1) = ',seriesname,';'];
    recurlist = [recurlist; {nextterm}];  
end

% Write MATLAB .m program file
funcname = func2str(fhandle);
if funcname(1) == '@' % anonymous function
    fileID1 = fopen('tempDEfunc.m','wt');
    fprintf(fileID1,'%s\n',['function dydt = tempDEfunc',funcname(2:6)]);
    fprintf(fileID1,'%s\n',['dydt = ',funcname(7:end),';']);
    fclose(fileID1);
    funcname = 'tempDEfunc';
end
newfuncname = [funcname,'series'];
newfunc = str2func(newfuncname);
fileID = fopen([newfuncname,'.m'],'wt');
headerline = ['function coefs = ',funcname,'series','(t0, y0, deg)'];
fprintf(fileID,'%s\n',headerline);
fprintf(fileID,'%s\n',['%',upper(newfuncname),' finds series coefs for soln of y'' = f(t,y) about t0 to deg']);

fprintf(fileID,'%s\n','% Where f(t,y) was given by:');
origfunctext = fileread([funcname,'.m']);
linecells = splitlines(origfunctext);
indent = '%    ';
newlinecells = cellfun(@(x)[indent,x],linecells,'UniformOutput',false);
fprintf(fileID,'%s\n',newlinecells{:});

fprintf(fileID,'%s\n',['% Automatically generated by makepsmcode.m, ',datestr(datetime),'.']);
%fprintf(fileID,'%% Equivalent polynomial DE system has max poly power m = %d; poly coef 1-norm C = %f.\n',maxm, maxC);

fopsdef = '  % number of scalar operations used to evaluate f';
fprintf(fileID,'\n%s\n',['numfops = ',num2str(length(valuelist)),';',fopsdef]);
fprintf(fileID,'%s\n','numys = length(y0);');
fprintf(fileID,'%s\n','T = [t0,1,zeros(1,deg-1)];');
fprintf(fileID,'%s\n','Y = zeros(numys,deg+1);');
fprintf(fileID,'%s\n','U = zeros(numfops,deg+1);');
fprintf(fileID,'%s\n','if isa(t0,''sym'') || isa(y0,''sym'')   % allow symbolic, including vpa');
fprintf(fileID,'%s\n','    T = sym(T); Y = sym(Y); U = sym(U);');
fprintf(fileID,'%s\n','end');
fprintf(fileID,'%s\n','Y(:,1) = y0;');

fprintf(fileID,'\n%s\n','% Evaluate f(t0,y0) recording value after each scalar operation');
fprintf(fileID,'%s\n',oplist{:});

fprintf(fileID,'\n%s\n','% Now recurrence rules for each operation in evaluation of f');
fprintf(fileID,'%s\n','for j = 2:deg   % j-1 is the power of term at index position j');
for line = 1:length(recurlist)
    fprintf(fileID,'%s\n',['    ',recurlist{line}]);
end
fprintf(fileID,'%s\n','end');

fprintf(fileID,'\n%s\n','coefs = Y;');
fprintf(fileID,'%s','end');
fclose(fileID);
message = [newfuncname, '.m  created or overwritten in current folder.'];
disp(message)
end