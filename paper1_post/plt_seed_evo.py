import uxarray as ux
import numpy as np
import holoviews as hv
import panel as pn
from holoviews import opts
import os
import glob
from paths import ARCHRT, CAMGR

DIRS = [ARCHRT[1], ARCHRT[3]]
FILS = [['0005-08-06-64800', '0005-08-07-00000', '0005-08-07-21600', '0005-08-07-43200'],
        ['0006-09-29-64800', '0006-09-30-00000', '0006-09-30-21600', '0006-09-30-43200']]

CTRS = [(98.79, 18.01), (22.73, 20.78)]
DEGW = 15

PLEVS = np.concatenate((np.arange(990, 1010, 5), np.arange(1010, 1026, 1)))
clim_min = float(PLEVS.min())
clim_max = float(PLEVS.max())

# Column headings and panel letters mapped linearly to the 4 columns
COL_TITLES = ["-6 hr", "0 hr", "+6 hr", "+12 hr"]
PANEL_LETTERS = [["(a)", "(b)", "(c)", "(d)"],  # Row 1
                 ["(e)", "(f)", "(g)", "(h)"]]  # Row 2

ROW_SUPERTITLES = ['UNSEED ' + FILS[0][1], 'SEED ' + FILS[1][1]]

# Fixed Customization Hook
def left_align_title(plot, element):
    try:
        # plot.state is the actual Bokeh Figure object
        if hasattr(plot.state, 'title') and plot.state.title is not None:
            plot.state.title.align = 'left'
    except Exception as e:
        print(f"Bypassing title alignment hook modification: {e}")

def main():
    pths = [[glob.glob(os.path.join(ar, 'atm/hist/', f'*h1i.{fi}*.nc')) for fi in FILS[ii]] for ii, ar in enumerate(DIRS)]

    print(pths)

    dss = [ux.open_mfdataset(CAMGR, [pt[0] for pt in pths[ii]]) for ii, ar in enumerate(DIRS)]

    subs = [ds['PSL'].subset.bounding_box((CTRS[ii][0] - DEGW, CTRS[ii][0] + DEGW), (CTRS[ii][1] - DEGW, CTRS[ii][1] + DEGW)) for ii, ds in enumerate(dss)]

    # 1. Keep this as a clear opts.Image type
    custom_opts = opts.Image(
        cmap='plasma', 
        clim=(clim_min, clim_max),
        color_levels=list(PLEVS),
        colorbar=True,
        frame_width=250,
        frame_height=250,
        data_aspect=1,
        aspect='equal',
        hooks=[left_align_title],
        # This will cleanly push the next column away without needing pn.Row(spacing=...)
        colorbar_opts={
            'padding': 10,       
            'margin': 35,        # Bumped to 35 to give plenty of clearance for pressure labels
            'label_standoff': 8  
        },
    )


    custom_opts = opts.Image(
        cmap='inferno', 
        clim=(clim_min, clim_max),
        color_levels=list(PLEVS),
        colorbar=True,
        
        # 1. Reset colorbar options to a simpler standard structure
        colorbar_opts={
            'padding': 5,
            'label_standoff': 8
        },
        
        # 2. Add an explicit margin to the plot object to prevent Panel truncation
        # Format: (left, top, right, bottom)
        # 45 pixels of strict empty breathing space to the right of every image frame!
        margin=(0, 0, 45, 0),
        
        frame_width=250,
        frame_height=250,
        data_aspect=1,
        hooks=[left_align_title]
    )

    grid_rows = []

    for ii, sb in enumerate(subs):
        row_plots = []
        for jj, tt in enumerate(sb['time'].values):
            full_title = f"{PANEL_LETTERS[ii][jj]}  {COL_TITLES[jj]}"
            
            # 2. Chain the .opts() calls to completely avoid the Tuple parameter error
            p = (sb.sel(time=tt) / 100).plot(rasterize=True).opts(
                custom_opts
            ).opts(
                title=full_title,
                title_format='{title}'
            )
            row_plots.append(p)
        
        # Assemble row plots, isolating axis sharing completely
        hv_row = hv.Layout(row_plots).cols(4).opts(shared_axes=False, axiswise=True)
        
        # Create a beautiful Markdown supertitle for this row segment
        super_title = pn.pane.Markdown(
            f"## **{ROW_SUPERTITLES[ii]}**", 
            styles={'margin': '25px 0 5px 10px', 'font-family': 'sans-serif'}
        )
        
        # Attach the label directly to the top of the row
        grid_rows.append(pn.Column(super_title, hv_row))

    # Combine everything vertically into a final Panel deployment 
    final_layout = pn.Column(*grid_rows)

    # Save to your interactive file
    final_layout.save('seed_evo_final.html')
    print("Interactive 2x4 grid with supertitles and customized headers saved successfully!")

if __name__ == "__main__":
    main()
