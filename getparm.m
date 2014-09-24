% function p=getparm(struct,parmname,parmdef);
function p=getparm(struct,parmname,parmdef);

if isfield(struct,parmname),
   p=getfield(struct,parmname);
elseif strcmp(parmdef,'ERROR'),
   error(sprintf('parameter %s required!',parmname));
   p=[];
else
   p=parmdef;
end
