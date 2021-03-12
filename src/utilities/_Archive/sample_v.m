function v_prop = sample_v(v_curr,sigma_v)
    v_prop = zeros(size(v_curr));
    for k = 1:numel(v_curr)
        pd = makedist('Normal',v_curr(k),sigma_v);
        t = truncate(pd,0,Inf);
        v_prop(k) = random(t);
    end
        