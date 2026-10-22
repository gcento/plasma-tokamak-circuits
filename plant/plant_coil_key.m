function key = plant_coil_key(coil)
% short coil key (e.g. 'CS3U', 'VS3U') of a pf_active.coil{i} element
% DD3: identifier when present and non-empty, else name
% DD4: name; a long name 'Long description (KEY)' gives 'KEY'
% Shared by configure_linear_model and the MUSCLE3 actor.

if isfield(coil,'identifier') && ~isempty(strtrim(coil.identifier))
    key = coil.identifier;
else
    key = coil.name;
    tok = regexp(key,'\(([^()]*)\)','tokens','once');
    if ~isempty(tok)
        key = tok{1};
    end
end
key = strtrim(key);
if isempty(key)
    error('linear_model:coilName','pf_active coil without name/identifier.');
end

end
