function signal = recover_signal(input, timeDelay, point, win_size, mean_power)

    tau = round(timeDelay / 2);
    win = floor(win_size / 2);

    a = max(1, point - win);
    b = min(length(input) - tau, point + win);

    idx = a:b;

    signal = input(idx + tau) ./ sqrt(mean_power);

end
