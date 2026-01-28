function power = meanPower(sig1, sig2, win_size, point)

    win = floor(win_size/2);

    a = max(1, point - win);
    b = min(length(sig1), point + win);

    N = b - a + 1;

    p1 = sum(sig1(a:b).^2) / N;
    p2 = sum(sig2(a:b).^2) / N;

    power = sqrt(p1 * p2);

end