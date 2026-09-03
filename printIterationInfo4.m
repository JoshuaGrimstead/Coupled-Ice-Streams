function printIterationInfo4(x)
    fprintf('Iter %d | T_s2 = %.2f km | W2 = %.3f\n', ...
        x.iter, x.T_s2, x.L1/1000);
end