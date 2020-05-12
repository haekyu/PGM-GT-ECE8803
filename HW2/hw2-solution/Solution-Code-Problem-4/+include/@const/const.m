classdef const < include.potential
% const Properties:
% value: numeric scalar
    methods
        function out=const(value)
            % pot=const(<value>)
            % value is a numeric scalar
            if nargin>0
                out.table=value;
                out.variables=[];
            end
        end
        
        %converters:
        
         function out=include.logconst(obj)  % convert const-> logconst
            out=include.logconst;
            out.table=log(obj.table);
        end
       
        function out=include.array(obj)  % convert const-> array
            out=include.array;
            out.variables=obj.variables;
            out.table=obj.table;
        end
        
        function out=include.logarray(obj)  % convert const-> logarray
            out=include.logarray;
            out.variables=obj.variables;
            out.table=log(obj.table);
        end
        
        function out=include.GaussianMoment(obj)  % convert const-> GaussianMoment
            out=include.GaussianMoment;        
            out.table.logprefactor=log(obj.table);
        end
        
    end
end