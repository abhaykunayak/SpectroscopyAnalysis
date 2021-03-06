classdef SpectroscopyData < handle
    properties
        hf1;
        p;
        X;
        Y;
        Z;
        V;
        LS;
        LS_com;
        LS_realspace;
        qx;
        qy
        I;
        I_cropped;
        I_slice;
        X_cropped;
        Y_cropped;
        Z_cropped;
        LS_cropped;
        LS_fft;
        LS_realspace_cropped;
        I_slice_cropped;
        ctr;
        hi1;
        ax1;
        hi2;
        ax2;
        hp_avg_ps;
        ax_avg_ps;
        hp_eline;
        slider;
        crop_button;
        crop_pos;
        reset_button;
        roi_imrect_handle;
        data;
        energy_q_axes_phase;
        energy_q_im_phase;
        LS_fft_phase;
    end
    methods
        
    end
end