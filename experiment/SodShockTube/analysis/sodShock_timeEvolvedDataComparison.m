function result = sodShock_timeEvolvedDataComparison(tubeLength)

    items       = evalin('base', 'who');
    mass        = zeros(tubeLength, 1);
    velocity    = zeros(tubeLength, 1);
    sorting     = zeros(length(items), 2);
    index       = 0;
    
    for i=1:length(items)
        if ~strcmpi(items{i}(1:3), 'sx_')
            fprintf('\nSkipped: %s', items{i});
            continue;
        end
        
        item                = evalin('base', items{i});
        index               = index + 1;
        time                = item.time.time;
        sorting(index, :)   = [index, time];
        s                   = SodShockSolution(length(item.mass), time);
        mass(:, index)      = abs(s.mass' - item.mass)./s.mass';
        velocity(:, index)  = abs(abs(s.velocity' - item.momX./item.mass)./s.velocity');
    end

    result.mass     = zeros(size(mass));
    result.velocity = zeros(size(velocity));
    sorting         = sortrows(sorting(1:index, :), 2);
    
    for i=1:index
        result.mass(:, i) = mass(:, sorting(i, 1));
        result.velocity(:, i) = velocity(:, sorting(i, 1));
    end
    
end