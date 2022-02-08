function filtered_proj=ram_lak_filter(projections,num_proj)
    %create a ram-lak filter
    N1 = size(projections,1);
    freqs=linspace(-1, 1, N1).';
    ram_filter = abs( freqs );
    ram_filter = repmat(ram_filter, [1 num_proj]);
    
    %filter each projection
    ft_R = fftshift(fft(projections,[],1),1);
    filteredProj = ft_R .* ram_filter;
    filteredProj = ifftshift(filteredProj,1);
    filtered_proj = real(ifft(filteredProj,[],1));
end
