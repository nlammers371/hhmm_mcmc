function fluo_out = fluo_conv_fun(initiation_rates,coeff_MS2)

    if size(initiation_rates,1) == 1
        initiation_rates = reshape(initiation_rates,[],1);
    end
    fluo_fragment = NaN(size(initiation_rates,1)+size(coeff_MS2,1)-1,size(initiation_rates,2),size(initiation_rates,3));
    for n = 1:size(initiation_rates,2)
        fluo_fragment(:,n,:) = convn(coeff_MS2(:,n),initiation_rates(:,n,:),'full');             
    end
    fluo_out = fluo_fragment(1:end-size(coeff_MS2,1)+1,:,:);