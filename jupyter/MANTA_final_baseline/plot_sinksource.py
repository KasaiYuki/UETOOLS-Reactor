from matplotlib.pyplot import subplots
from matplotlib.figure import Figure


def plot_sinksource(case):
    f = case.plot.grid(plates=False, vessel=False)
    plot_sinks(case, ax=f)
    plot_sources(case, ax=f)


def plot_sinks(case, ax=None, maxpump=0.3):
    if ax is None:
        f, ax = subplots()
    elif isinstance(ax, Figure):
        f = ax
        ax = f.get_axes()[0]
    else:
        f = ax.get_figure()
    colors = ['b', 'c', 'lightblue', 'teal']

    rm = case.get('rm')
    zm = case.get('zm')
    albi = case.get('albedoi')
    albo = case.get('albedoo')
    albl = case.get('alblb')
    albr = case.get('albrb')
    ixrb = case.get('ixrb')
    ixlb = case.get('ixlb')
    # Loop over species
    for s in range(case.get('ngsp')):
        # Inner radial boundary
        for ix in range(len(albi)):
            ax.plot(rm[ix, 1, 1:3], zm[ix, 1, 1:3], color=colors[s],
                linewidth=3, alpha = min((1-albi[ix,s])/maxpump, 1))
        # Outer radial boundary
        for ix in range(len(albo)):
            ax.plot(rm[ix, -1, 1:3], zm[ix, -1, 1:3], color=colors[s],
                linewidth=3, alpha = min((1-albo[ix,s])/maxpump, 1))
        for p in range(len(ixlb)):
            for iy in range(len(albl)):
                ax.plot(rm[ixlb[p], iy, 2:4], zm[ixlb[p], iy, 2:4], 
                    color=colors[s], linewidth=3, 
                    alpha = min((1-albl[iy,s,p])/maxpump, 1)
                )
            for iy in range(len(albr)):
                ax.plot(rm[ixrb[p], iy, 2:4], zm[ixrb[p], iy, 2:4], 
                    color=colors[s], linewidth=3, 
                    alpha = min((1-albr[iy,s,p])/maxpump, 1)
                )

def plot_sources(case, ax=None):
    if ax is None:
        f, ax = subplots()
    elif isinstance(ax, Figure):
        f = ax
        ax = f.get_axes()[0]
    else:
        f = ax.get_figure()
    colors = ['r', 'm', 'orange', 'gold']

    rm = case.get('rm')
    zm = case.get('zm')
    fngyso = case.get('fngyso')
    fngysi = case.get('fngysi')
    maxpuff = max(fngyso.max(), fngysi.max())
    # Loop over species
    for s in range(case.get('ngsp')):
        # Inner radial boundary
        for ix in range(len(fngysi)):
            ax.plot(rm[ix, 1, 1:3], zm[ix, 1, 1:3], color=colors[s],
                linewidth=3, alpha = fngysi[ix,s]/maxpuff)
        # Outer radial boundary
        for ix in range(len(fngyso)):
            ax.plot(rm[ix, -1, 1:3], zm[ix, -1, 1:3], color=colors[s],
                linewidth=3, alpha = fngyso[ix,s]/maxpuff)
