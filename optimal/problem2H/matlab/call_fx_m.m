function [fx gr] = call_fx_m(x)
    % Calculate objective fx
    % disp('Calculating function value.');
    fx = calllib('problem2H', 'call_fx', x);
    % disp('Function calculated; '+fx);
    % fprintf('Function calculated %f\n',fx);
    
    if nargout > 1 %gradient required
        % disp('---> Calculating gradient value.');
        gr = zeros(1,16);
        [x gr] = calllib('problem2H', 'call_gr', x, gr, 16);
        % disp(x);
        % disp(gr);
        % fprintf('rx: '); fprintf('%8.4f ', x); fprintf('\n');
        % fprintf('gr: '); fprintf('%8.4f ', gr); fprintf('\n');
        % disp('---> Calculate gradients');
    end
    
end