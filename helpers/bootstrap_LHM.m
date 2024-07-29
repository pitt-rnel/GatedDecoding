function[LHM] = bootstrap_LHM(func,X)
    
    [bl,bh] = boot_bounds(1000,func,X,2.5,97.5); 
    
    bm = func(X); 
    if size(bm,2)>size(bm,1)
        bm = bm(:); 
    end

    LHM = [bl,bh,bm]; 
end