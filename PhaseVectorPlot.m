function [] = PhaseVectorPlot()
    ax = gca;
    im = imhandles(ax);
    
    figure;
    ax_phase = axes;
    q_phase = quiver([],[],[],[]);
    
    roi = drawrectangle(ax);
    addlistener(roi,'MovingROI',@(src, evt) updatePhasePlot(src, evt,...
        q_phase, ax_phase, im, ax));
    
end

function updatePhasePlot(src, evt, q_phase, ax_phase, im, ax)
    mag_data = ax.UserData;
    data = im.CData;
    xdata = im.XData;
    ydata = im.YData;
    [nrows, ncols] = size(data);
    
    pos = evt.CurrentPosition;
    px = int16(axes2pix(ncols, xdata, [pos(1) pos(1)+pos(3)]));
    py = int16(axes2pix(nrows, ydata, [pos(2) pos(2)+pos(4)]));
    
    data_crop = data(py(2):py(1), px(1):px(2));
    xdata_crop = xdata(px(1):px(2));
    ydata_crop = ydata(py(2):py(1))*100;
    r = mag_data(py(2):py(1), px(1):px(2));
%     r = 1;
    set(q_phase, 'XData', xdata_crop);
    set(q_phase, 'YData', ydata_crop);
    set(q_phase, 'UData', r.*cos(data_crop));
    set(q_phase, 'VData', r.*sin(data_crop));
    axis(ax_phase, 'tight');
end
    
