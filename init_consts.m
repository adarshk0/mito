    function ret = init_consts(input, mito_init, myo_init, im)
        ret = input;
        tmp = ret(:,:,1);
        tmp(im)  =  mito_init;
        tmp(~im) =  myo_init;
        ret(:,:,1) = tmp;
    end