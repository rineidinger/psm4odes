classdef trace
   %trace class: literal & numerical lists of expression evaluation steps.
   % Designed for expressions of form f(t,y) (for use in y' = f(t,y)) where
   %    t is independent variable, construct: t = trace(t0, 0)
   %    y is a vector variable (nx1), construct: y = trace(Y0, 1)
   % Every overload operation (standard scalar functions and operations)
   % on a trace object keeps operation count and appends to lists or maxs:
   %    count: integer (in double) for current operation count
   %    vallist: column matrix (countx1) of all computation step doubles
   %    oplist: cell array of strings defining each intermediate
   %    reclist: cell array of strings relating how series can be combined
   % After evaluation, call
   %    [count, vallist, oplist, reclist] = trace.count_val
   % by Richard Neidinger 1/11/21; rev. 10/26/22, formatting 7/1/24,
   % finalized for publication 7/19/24.
   properties
      % One trace object simply holds the current value and operation count
      val  % expression value after its operation/function
      opcount % index of val in vallist maintained by count_val
      tlinear % T/F: operations to this point are equivalent to at+b
              %   Boolean true=1 if series coefs of order >1 are zero. Will
              %   simplify code where dot product has only one nonzero term
   end 
   
   methods (Static)
       function [op, val, code, rec]...
               = count_val(value, opcell, reccell)
           % trace.count_val will append and/or return count and lists
           persistent count vallist oplist reclist
           if isempty(count)  % on first call, all are initialized
               % default intialization works for vallist = [];
               count = 0;
               oplist = {};
               reclist = {};
           end
           if nargin % any argument, so assignment, not just returns
               count = count + 1;
               vallist = [vallist; value];
               if nargin > 1
                   oplist = [oplist; opcell];
                   reclist = [reclist; reccell];
               end
           end
           op = count;  % return count and all persistent variable values
           val = vallist;
           code = oplist;
           rec = reclist;
       end
       function opstr = name(i)
           % returns string 'T(1', 'U(5', or 'Y(5' for i = 0, 5, -5 resp.
           if i == 0
               opstr = 'T(1';
           elseif i < 0
               opstr = ['Y(',num2str(-i)];
           else
               opstr = ['U(',num2str(i)];
           end
       end
       function opstr = nextname(advance)
           % returns string ui where i is opcount + advance
           % advance is an optional argument, assumed 1
           if ~nargin
               advance = 1;
           end
           opstr = ['U(',num2str(advance+trace.count_val)];
       end
   end
   
   methods
      function obj = trace(value, role)
         %trace class constructor for value of variable or operation
         % role: 0 is independent variable (called t)
         %       1 is dependent DE vector variable (called y)
         %       2 is indermediate variable (called u with opcount index)
         if nargin == 0  % Create empty fields to be assigned before use.
             obj.val = [];
             obj.opcount = NaN; % see next case 
             obj.tlinear = [];  
         elseif nargin == 1 || role == 2  % intermediate variable default
             obj.val = value;
             obj.opcount = NaN; % set by count_val in calling function
             % with no calling operation, NaN flags a pure constant
             obj.tlinear = []; % set by calling function
         elseif role == 0
             obj.val = value;
             obj.opcount = 0; % indicates independent variable
             obj.tlinear = true;
         elseif role == 1
             for i = 1:length(value)  % so each element of y vector is obj
                % obj(i) implicitly calls the nargin == 0 case above
                obj(i).val = value(i);
                obj(i).opcount = -i; % negatives indicate "y" variable
                obj(i).tlinear = false;
             end
         else % not intended for use
             disp('Warning: trace called with role other than 0, 1, 2')
             obj.val = [];
             obj.opcount = [];
             obj.tlinear = [];
             obj.m = [];
             obj.C = [];
         end
      end
      function [val, opcount] = valueandopcount(obj)
         %trace/VALUEANDOPCOUNT Vector value(s) and opcount(s) of object(s)
         % was DOUBLE method but renamed so double = object won't convert
         val = [obj.val].';
         opcount = [obj.opcount].'; 
      end
      function h = plus(u,v)
         %trace/PLUS overloads addition with at least one trace arg
         % h = u + v  for both operation value and series recurrence code
         if ~isa(u,'trace') % assume u is a numerical scalar
             newval = u + v.val;
             vname = trace.name(v.opcount);
             opstr = [num2str(u,16),' + ',vname,',1);'];
             recstr = [vname,',j);'];
             newtlinear = v.tlinear;
         elseif ~isa(v,'trace') % assume v is a numerical scalar
             newval = u.val + v;
             uname = trace.name(u.opcount);
             opstr = [uname,',1) + ',num2str(v,16),';'];
             recstr = [uname,',j);'];
             newtlinear = u.tlinear;
         else  
             newval = u.val + v.val;
             uname = trace.name(u.opcount);
             vname = trace.name(v.opcount);
             opstr = [uname,',1) + ',vname,',1);'];
             recstr = [uname,',j) + ',vname,',j);'];
             newtlinear = v.tlinear & u.tlinear;
         end
         h = trace(newval);
         hname = trace.nextname;
         opstr  = [hname,',1) = ',opstr];
         recstr = [hname,',j) = ',recstr];
         h.opcount = trace.count_val(newval, {opstr}, {recstr});
         h.tlinear = newtlinear; 
      end
      function h = uminus(u)
         %trace/UMINUS overloads negation with a trace object argument
         % h = -u  for both operation value and series recurrence code
         newval = -u.val;
         h = trace(newval);
         opstr  = [trace.nextname,',1) = -', trace.name(u.opcount),',1);'];
         recstr = [trace.nextname,',j) = -', trace.name(u.opcount),',j);'];
         h.opcount = trace.count_val(newval, {opstr}, {recstr});
         h.tlinear = u.tlinear;
      end
      function h = minus(u,v)
         %trace/MINUS overloads subtraction with trace argument(s)
         % h = u - v  for both operation value and series recurrence code
         if ~isa(u,'trace') % assume u is a numerical scalar
             newval = u - v.val;
             vname = trace.name(v.opcount);
             opstr = [num2str(u,16),' - ',vname,',1);'];
             recstr = ['-',vname,',j);'];
             newtlinear = v.tlinear;
         elseif ~isa(v,'trace') % assume v is a numerical scalar
             newval = u.val - v;
             uname = trace.name(u.opcount);
             opstr = [uname,',1) - ',num2str(v,16),';'];
             recstr = [uname,',j);'];
             newtlinear = u.tlinear;
         else  
             newval = u.val - v.val;
             uname = trace.name(u.opcount);
             vname = trace.name(v.opcount);
             opstr = [uname,',1) - ',vname,',1);'];
             recstr = [uname,',j) - ',vname,',j);'];
             newtlinear = v.tlinear & u.tlinear;
         end
         opstr  = [trace.nextname,',1) = ',opstr];
         recstr = [trace.nextname,',j) = ',recstr];
         h = trace(newval);
         h.opcount = trace.count_val(newval, {opstr}, {recstr});
         h.tlinear = newtlinear;
      end
      function h = mtimes(u,v)
         %trace/MTIMES overloads multiplication with trace arg(s)
         % h = u * v  gives value and Cauchy product for recurrence code
         hname = trace.nextname;
         if ~isa(u,'trace') % assume u is a scalar
             newval = u * v.val;
             midstr = [num2str(u,16),' * ',trace.name(v.opcount)];
             opstr  = [hname,',1) = ',midstr,',1);'];
             recstr = [hname,',j) = ',midstr,',j);'];
             newtlinear = v.tlinear;
         elseif ~isa(v,'trace')
             newval = u.val * v;
             midstr = [num2str(v,16),' * ',trace.name(u.opcount)];
             opstr =  [hname,',1) = ',midstr,',1);'];
             recstr = [hname,',j) = ',midstr,',j);'];
             newtlinear = u.tlinear;
         else
             newval = u.val * v.val;
             uname = trace.name(u.opcount);
             vname = trace.name(v.opcount);
             opstr = [hname,',1) = ',uname,',1) * ',vname,',1);'];
             recstr = [hname,',j) = ',uname,',j:-1:1) * ',vname,',1:j).'';'];
             newtlinear = false;
         end
         h = trace(newval, 2);
         h.opcount = trace.count_val(newval, {opstr}, {recstr});
         h.tlinear = newtlinear;  
      end
      function h = mrdivide(u,v)
         %trace/MRDIVIDE overloads division with trace arg(s)
         % h = u/v gives value; for recurrence code:
         % u = h * v  gives Cauchy product solved for highest h coef
         hname = trace.nextname;
         if ~isa(u,'trace') % assume u is a scalar
             newval = u / v.val;
             vname = trace.name(v.opcount);
             opstr  = [hname,',1) = ',num2str(u,16),' / ',vname,',1);'];
             if v.tlinear
                 recstr = [hname,',j) = -',hname,',j-1) * ',vname,',2) / ',vname,',1);'];
             else
                 recstr = [hname,',j) = -( ',vname,',j:-1:2) * ',hname,',1:(j-1)).'' ) / ',vname,',1);'];
             end
             newtlinear = false;
         elseif ~isa(v,'trace')
             newval = u.val / v;
             uname = trace.name(u.opcount);
             opstr =  [hname,',1) = ',uname,',1) / ',num2str(v,16),';'];
             recstr = [hname,',j) = ',uname,',j) / ',num2str(v,16),';'];
             newtlinear = u.tlinear;
         else
             newval = u.val / v.val;
             uname = trace.name(u.opcount);
             vname = trace.name(v.opcount);
             opstr = [hname,',1) = ',uname,',1) / ',vname,',1);'];
             recstart = [hname,',j) = ( ',uname,',j) - '];
             if v.tlinear
                 recstr = [recstart,hname,',j-1) * ',vname,',2) ) / ',vname,',1);'];
             else
                 recstr = [recstart,vname,',j:-1:2) * ',hname,',1:(j-1)).'' ) / ',vname,',1);'];
             end
             newtlinear = false;
         end
         h = trace(newval, 2);
         h.opcount = trace.count_val(newval, {opstr}, {recstr});
         h.tlinear = newtlinear;
      end
      function h = sqrt(u)
         %trace/SQRT overloads square root of a trace object argument
         % h = sqrt(u) gives value; for recurrence code:
         % u = h * h  gives Cauchy product solved for highest h coef
         newval = sqrt(u.val);
         h = trace(newval, 2);
         hname = trace.nextname;
         uname = trace.name(u.opcount);
         opstr = [hname,',1) = sqrt(', uname,',1));'];
         recstr = [hname,',j) = ( ',uname,',j) - ',hname,',(j-1):-1:2) * ',hname,',2:(j-1)).'' ) / (2*',hname,',1));'];
         h.opcount = trace.count_val(newval, {opstr}, {recstr});
         h.tlinear = false;
      end
      function h = exp(u)
         %trace/EXP overloads exp with a trace object argument
         % h = exp(u) gives value; for recurrence code:
         % h' = h * u'  gives Cauchy product solved for highest h coef
         newval = exp(u.val);
         h = trace(newval, 2);
         hname = trace.nextname;
         uname = trace.name(u.opcount);
         opstr = [hname, ',1) = exp(', uname,',1));'];         
         if u.tlinear
            reccell = {[hname,',j) = ',hname,',j-1) * ',uname,',2) /(j-1);']};
         else
            % tempprime pulls coefs for u' out of array coefs for u
            recstr1 = ['tempprime = ( ',uname,',2:j) .* (1:(j-1)) ).'';'];
            recstr2 = [hname,',j) = ( ',hname,',(j-1):-1:1) * tempprime )/(j-1);'];
            reccell = {recstr1;recstr2};
         end
         h.opcount = trace.count_val(newval, {opstr}, reccell);
         h.tlinear = false;
      end
      function h = log(u)
         %trace/LOG overloads natural logarithm of a trace argument
         % h = log(u) gives value; for recurrence code:
         % u' = u * h'  gives Cauchy product solved for highest h coef
         newval = log(u.val);
         h = trace(newval, 2);
         hname = trace.nextname;
         uname = trace.name(u.opcount);
         opstr = [hname, ',1) = log(', uname,',1));'];
         if u.tlinear
            recpart = [uname,',2)*',hname,',j-1)/(j-1)'];
            recstr = [hname,',j) = ( ',uname,',j) - ( (j-2)*',recpart,' ) ) /',uname,',1);'];
            reccell = {recstr};
         else
            recstr1 = ['tempprime = ( ',hname,',2:(j-1)) .* (1:(j-2)) ).'';'];
            recstr2 = [hname,',j) = ( ',uname,',j) - (',uname,',(j-1):-1:2) * tempprime )/(j-1) ) /',uname,',1);'];
            % In case j=2, vectors are empty and * results in 0
            reccell = {recstr1;recstr2};
         end
         h.opcount = trace.count_val(newval, {opstr}, reccell);
         h.tlinear = false;
      end
      function h = mpower(u,r)
         %trace/MPOWER overloads power with at least one trace arg
         if isa(r,'trace')
             h = exp(r*log(u));  % if u is scalar, log(u) is just a scalar
         else % assume r is a scalar
             if r == 2 % first special cases use only one Cauchy product
                 h = u*u;
             elseif r == 1/2
                 h = sqrt(u);
             elseif r == -1
                 h = 1/u;
             elseif r == 3 % two Cauchy products but no domain restriction
                 h = u*u*u;
             elseif r == 4
                 v = u*u;
                 h = v*v;
             else % general r scalar allows nonzero neg base (unlike exp&log)
                 % h = u^r gives value and h' = r*u^(r-1)*u', so 
                 % u*h' = r*h*u' gives Cauchy products solved for highest h coef
                 newval = u.val ^ r;
                 h = trace(newval, 2);
                 hname = trace.nextname;
                 uname = trace.name(u.opcount);
                 opstr = [hname,',1) = ',uname,',1)^',num2str(r,16),';'];
                 if u.tlinear
                    recstart = [hname,',j) = (',num2str(r,16),'-(j-2))*'];
                    reccell = {[recstart,uname,',2)*',hname,',j-1) / ((j-1)*',uname,',1));']};
                 else
                    uprimestr = ['tempinprime = ( ',uname,',2:j) .* (1:(j-1)) ).'';'];
                    hprimestr = ['tempoutprime = ( ',hname,',2:(j-1)) .* (1:(j-2)) ).'';'];
                    cp1str = ['cp1 = ',hname,',(j-1):-1:1) * tempinprime;'];
                    cp2str = ['cp2 = ',uname,',(j-1):-1:2) * tempoutprime;'];
                    recstr = [hname,',j) = ( ',num2str(r,16),'*cp1  - cp2 ) / ( (j-1)*',uname,',1) );'];
                    reccell = {uprimestr; hprimestr; cp1str; cp2str; recstr};
                 end
                 h.opcount = trace.count_val(newval, {opstr}, reccell);
                 h.tlinear = false;
             end
         end
      end
      function [s, c] = sincos(u)
         %trace/SINCOS returns trace objects for sin(u) and cos(u)
         % s = sin(u) and c = cos(u) 
         % s' = c * u' and c' = -s * u' Cauchy products series recurrence
         % Code for above four caluculations are written into strings using
         % next array&index names for s and c and array&index name for u.
         % Strings are put into cells and appended to oplist and reclist by
         % trace.count_val which also updates operation count.
         newval = sin(u.val);
         s = trace(newval);
         sname = trace.nextname;    % 'U(i' where i is trace.count_val + 1
         cname = trace.nextname(2); % 'U(i' where i is trace.count_val + 2
         uname = trace.name(u.opcount);
         opstr = [sname, ',1) = sin(', uname,',1));'];
         if u.tlinear
            reccell = {[sname,',j) = ',cname,',j-1)*',uname,',2)/(j-1);']};
         else
            % tempprime pulls coefs for u' out of array coefs for u
            recstr1 = ['tempprime = ( ',uname,',2:j) .* (1:(j-1)) ).'';'];
            recstr2=[sname,',j) = (',cname,',(j-1):-1:1) * tempprime)/(j-1);'];
            reccell = {recstr1;recstr2};
         end
         s.opcount = trace.count_val(newval, {opstr}, reccell);
         s.tlinear = false;
         
         newval = cos(u.val);
         c = trace(newval);
         opstr = [cname, ',1) = cos(', uname,',1));'];
         if u.tlinear
            recstr = [cname,',j) = -',sname,',j-1)*',uname,',2)/(j-1);'];
         else
            recstr=[cname,',j) = -(',sname,',(j-1):-1:1) * tempprime)/(j-1);'];
         end
         c.opcount = trace.count_val(newval, {opstr}, {recstr});
         c.tlinear = false;    
      end
      function h = sin(u)
         %trace/SIN overloads sine with a trace object argument
         [h, g] = sincos(u);
      end
      function g = cos(u)
         %trace/COS overloads cosine of a trace object argument
         [h, g] = sincos(u);
      end
      function h = tan(u)
         %trace/TAN overloads tangent of a trace object argument
         % h = tan(u) and v = sec^2(u) = 1 + h^2 give values
         % h' = v*u' gives Cauchy product solved for h recurrence
         % v' = 2*h*h' gives Cauchy product solved for v recurrence
         newval = tan(u.val);
         h = trace(newval, 2);
         hname = trace.nextname;
         vname = trace.nextname(2);
         uname = trace.name(u.opcount);
         opstr = [hname,',1) = tan(', uname,',1));'];
         if u.tlinear
            reccell = {[hname,',j) = ',vname,',j-1) * ',uname,',2) /(j-1);']};
         else
            % tempprime pulls coefs for u' out of array coefs for u
            recstr1 = ['tempprime = ( ',uname,',2:j) .* (1:(j-1)) ).'';'];
            recstr2 = [hname,',j) = (',vname,',(j-1):-1:1) * tempprime)/(j-1);'];
            reccell = {recstr1;recstr2};
         end
         h.opcount = trace.count_val(newval, {opstr}, reccell);
         h.tlinear = false;

         newval = 1 + u.val^2;
         v = trace(newval, 2);  % object probably never used, just atrings
         opstr = [vname, ',1) = 1 + ', hname,',1)^2;'];
         recstr1 = ['tempprime = ( ',hname,',2:j) .* (1:(j-1)) ).'';'];
         recstr2 = [vname,',j) = 2*( ',hname,',(j-1):-1:1) * tempprime )/(j-1);'];
         reccell = {recstr1;recstr2};
         v.opcount = trace.count_val(newval, {opstr}, reccell);
         v.tlinear = false;   
      end
      function h = cot(u)
         %trace/COT overloads tangent of a trace object argument
         % h = tan(u) and v = csc^2(u) = 1 + h^2 give values
         % h' = -v*u' gives Cauchy product solved for h recurrence
         % v' = 2*h*h' gives Cauchy product solved for v recurrence
         newval = tan(u.val);
         h = trace(newval, 2);
         hname = trace.nextname;
         vname = trace.nextname(2);
         uname = trace.name(u.opcount);
         opstr = [hname,',1) = cot(', uname,',1));'];
         if u.tlinear
            reccell = {[hname,',j) = - ',vname,',j-1) * ',uname,',2) /(j-1);']};
         else
            % tempprime pulls coefs for u' out of array coefs for u
            recstr1 = ['tempprime = ( ',uname,',2:j) .* (1:(j-1)) ).'';'];
            recstr2 = [hname,',j) = - (',vname,',(j-1):-1:1) * tempprime)/(j-1);'];
            reccell = {recstr1;recstr2};
         end
         h.opcount = trace.count_val(newval, {opstr}, reccell);
         h.tlinear = false;

         newval = 1 + u.val^2;
         v = trace(newval, 2);  % object probably never used, just atrings
         opstr = [vname, ',1) = 1 + ', hname,',1)^2;'];
         recstr1 = ['tempprime = ( ',hname,',2:j) .* (1:(j-1)) ).'';'];
         recstr2 = [vname,',j) = 2*( ',hname,',(j-1):-1:1) * tempprime )/(j-1);'];
         reccell = {recstr1;recstr2};
         v.opcount = trace.count_val(newval, {opstr}, reccell);
         v.tlinear = false;   
      end
      function h = asin(u)
         %trace/ASIN overloads arcsine of a trace object argument
         % h = asin(u) and v = cos(h) give values,  since u = sin(h),
         % u' = v*h' gives Cauchy product solved for h recurrence
         % v' = -u*h' gives Cauchy product recurrence for v
         newval = asin(u.val);
         h = trace(newval, 2);
         hname = trace.nextname;
         vname = trace.nextname(2);
         uname = trace.name(u.opcount);
         opstr = [hname, ',1) = asin(', uname,',1));'];
         recstr1 = ['tempprime = ( ',hname,',2:(j-1)) .* (1:(j-2)) ).'';'];
         recstr2 = [hname,',j) = ( ',uname,',j) - (',vname,',(j-1):-1:2) * tempprime )/(j-1) ) /',vname,',1);'];
         % In case j=2, vectors are empty and * results in 0
         reccell = {recstr1;recstr2};
         h.opcount = trace.count_val(newval, {opstr}, reccell);
         h.tlinear = false;

         newval = cos(h.val);
         v = trace(newval, 2);
         opstr = [vname, ',1) = cos(', hname,',1));'];
         % not simplified for tlinear, different j=2 case handled by empty tempprime
         recstr = [vname,',j) = -',uname,',1)*',hname,',j) - (',uname,',(j-1):-1:2) * tempprime )/(j-1);'];
         v.opcount = trace.count_val(newval, {opstr}, {recstr});
         v.tlinear = false;
      end
      function h = atan(u)
         %trace/ATAN overloads arctangent of a trace object argument
         % h = atan(u) and v = 1 + u^2 = sec^2(h) give values
         % then u = tan(h) and u' = v*h'
         % u' = v*h' gives Cauchy product solved for h recurrence
         % v' = 2u*u' gives Cauchy product recurrence for v
         newval = atan(u.val);
         h = trace(newval, 2);
         hname = trace.nextname;
         vname = trace.nextname(2);
         uname = trace.name(u.opcount);
         opstr = [hname, ',1) = atan(', uname,',1));'];
         recstr1 = ['tempprime = ( ',hname,',2:(j-1)) .* (1:(j-2)) ).'';'];
         recstr2 = [hname,',j) = ( ',uname,',j) - (',vname,',(j-1):-1:2) * tempprime )/(j-1) ) /',vname,',1);'];
         % In case j=2, vectors are empty and * results in 0
         reccell = {recstr1;recstr2};
         h.opcount = trace.count_val(newval, {opstr}, reccell);
         h.tlinear = false;

         newval = 1 + u.val^2;
         v = trace(newval, 2);  
         opstr = [vname, ',1) = 1 + ', uname,',1)^2;'];
         if u.tlinear  % (only j=2 & j=3 nonzero)
            reccell = {[vname,',j) = 2 * ',uname,',(j-1)) * ',uname,',2) / (j-1);']};
         else
            recstr1 = ['tempprime = ( ',uname,',2:j) .* (1:(j-1)) ).'';'];
            recstr2 = [vname,',j) = 2 * (',uname,',(j-1):-1:1) * tempprime )/(j-1);'];
            reccell = {recstr1;recstr2};
         end
         v.opcount = trace.count_val(newval, {opstr}, reccell);
         v.tlinear = false;
      end
      % Many other functions can be combinations of above.
      function h = sec(u)
         %trace/SEC overloads sec of a trace object argument
         h = 1/cos(u);
      end 
      function h = csc(u)
         %trace/CSC overloads csc of a trace object argument
         h = 1/sin(u);
      end      
      function h = acos(u)
         %trace/ACOS overloads acos of a trace object argument
         h = pi/2 - asin(u);
      end       
      function h = sinh(u)
         %trace/SINH overloads sinh of a trace object argument
         h = (exp(u) - exp(-u))/2;
      end  
      function h = cosh(u)
         %trace/COSH overloads cosh of a trace object argument
         h = (exp(u) + exp(-u))/2;
      end
      % tanh and coth can be done directy using 2 Cauchy products
      % h = tanh(u) or h = coth(u), v = 1-h^2, so h' = v*u', v' = -2h*h'
      function h = atanh(u)
         %trace/ATANH overloads arctangent of a trace object argument
         h = (log(1+u) - log(1-u))/2;
      end
      % asinh and acosh can be done directy using 2 Cauchy products
      % h = asinh(u), v = cosh(h), so h' = u'/v, v' = u*h'
      % h = acosh(u), v = sinh(h), so h' = u'/v, v' = u*h'
   end
       
end