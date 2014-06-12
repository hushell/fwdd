function answer = isoctave()
% Return true if this function is run on Octave, otherwise false.     
    answer = ~isSubstring('MATLAB', matlabroot, true); 
end