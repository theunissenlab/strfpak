function varargout = whos_test(n)
for jj = 1:nargout
    eval(['out' num2str(jj) ' = rand(' num2str(jj) ');']);
end
trouble = rand(400);
whos_out = whos('out*');
sum([whos_out(:).bytes])
for jj = 1:nargout
    eval(['varargout{' num2str(jj) '} = out' num2str(jj) ';']);
end

save('whos_test_out','out*');

whos_2 = whos('varargout');
whos_2.bytes
