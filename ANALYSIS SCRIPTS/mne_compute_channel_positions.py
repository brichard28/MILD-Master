def mne_compute_channel_positions(info, n_channels):
    import numpy as np
    import mne
    from mne.viz.utils import (apply_trans, _check_sphere)
    from mne.viz.topomap import _get_pos_outlines
    from mne._fiff.constants import FIFF

    dev_head_t = info["dev_head_t"]
    chs = info["chs"][:]
    pos = np.empty((len(chs), 3))
    for ci, ch in enumerate(chs):
        pos[ci] = ch["loc"][:3]
        if ch["coord_frame"] == FIFF.FIFFV_COORD_DEVICE:
            if dev_head_t is None:
                dev_head_t = np.eye(4)
            pos[ci] = apply_trans(dev_head_t, pos[ci])
    del dev_head_t

    sphere = None
    sphere = _check_sphere(sphere, info)
    picks = np.arange(n_channels)
    pos, outlines = _get_pos_outlines(info, picks, sphere, to_sphere=True)

    return pos, outlines