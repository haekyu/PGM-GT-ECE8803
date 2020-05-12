function newpot=orderpot(pot,varargin)
%ORDERPOT Return potential with variables reordered according to order
% newpot=orderpot(pot,<order>)
% if order is missing or empty, the variables are sorted (low to high)
if nargin==2; inorder=varargin{1}; else inorder=[]; end
pot=include.str2cell(pot);

if length(pot)>1
    for p=1:length(pot)
        if isempty(inorder)
            order=sort(pot{p}.variables);
        else
            order=inorder;
        end
        
        [a iperm]=ismember(order,pot{p}.variables);
        iperm=iperm(iperm>0);
        
        if isa(pot{p},'include.const') || isa(pot{p},'include.logconst')
            newpot{p}=pot{p};
        else
            newpot{p}=permute(pot{p},iperm);
        end
    end
else
    if isempty(inorder)
        order=sort(pot{1}.variables);
    else
        order=inorder;
    end
    
    [a iperm]=ismember(order,pot{1}.variables);
    iperm=iperm(iperm>0);
    
    if isa(pot{1},'include.const') || isa(pot{1},'include.logconst')
        newpot=pot{1};
    else
        newpot=permute(pot{1},iperm);
    end
end