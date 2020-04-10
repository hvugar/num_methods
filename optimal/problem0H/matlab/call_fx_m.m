function [fx gr] = call_fx_m(x)
    % Calculate objective fx
%     disp('Calculating function value...');
%     disp(x);
    fx = calllib('problem0H', 'call_fx', x);
%     disp('Function calculated; ');
%     fprintf('Function calculated %f\n',fx);
    
    if nargout > 1 %gradient required
%         disp('---> Calculating gradient value.');
        gr = zeros(1,16);
        [x gr] = calllib('problem0H', 'call_gr', x, gr, 16);
        % disp(x);
        % disp(gr);
        % fprintf('rx: '); fprintf('%8.4f ', x); fprintf('\n');
        % fprintf('gr: '); fprintf('%8.4f ', gr); fprintf('\n');
        % disp('---> Calculate gradients');
    end
    
end