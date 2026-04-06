function signal = recover_signal(input, timeDelay, point, win_size)

    tau = round(timeDelay / 2); %It's divided 2 because we increase one and decrase the other (it's DToD)
    win = floor(win_size / 2);

    a = max(1, point - win);
    b = min(length(input) - tau, point + win);

    idx = a:b;

    signal = input(idx + tau); %We apply to all points in the window the same delay

end
