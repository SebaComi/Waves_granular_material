function save_GIF(salva, k, filename, Deltat)
if salva
    % Capture the frame
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to the GIF file
    if k == 1
        imwrite(imind, cm, ['Figure/' filename], 'gif', 'Loopcount', 0);%, 'DelayTime', 30*Deltat);
    else
        imwrite(imind, cm, ['Figure/' filename], 'gif', 'WriteMode', 'append');%, 'DelayTime', 30*Deltat);
    end
end