function out = get_ntrials(DS);
out = [];

for jj = 1:length(DS)
    out(jj) = DS{jj}.ntrials;
end
