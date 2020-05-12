classdef (InferiorClasses = {?include.const,?include.logconst}) array < include.potential
% array Properties:
% variables
% table
%
% Inferior classes: include.const, include.logconst
% type help array.array for constructor help    
    methods
        
        function arr=array(variables,table)
            % array(<variables>,<table>)
            % eg pot=array([1 2],rand(2,2))
            % or pot=array; pot.variables=[1 2]; pot.table=rand(2,2);
            if nargin>0
                if ~isnumeric(variables)
                    error('variables must be a numerical vector')
                else
                    arr.variables=variables;
                end
                
                if nargin>1
                    if ~isnumeric(table)
                        error('table must be a numerical array')
                    else
                        arr.table=table;
                    end
                    
                    if length(variables)~=length(include.mysize(table))
                        error('number of declared variables is not equal to the number of variables in the table')
                    end
                end
            end
        end
        
        % place converters here:
        function out=include.logarray(obj) % convert array -> logarray
            out=include.logarray;
            out.variables=obj.variables;
            out.table=log(obj.table);
        end
    end
end