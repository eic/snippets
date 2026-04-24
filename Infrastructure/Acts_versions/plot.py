from datetime import datetime, timedelta

import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

COMPONENT = "DD4hep"

@np.vectorize
def to_datetime(date):
    return datetime.strptime(date, '%b %d %Y %H:%M:%S') 

with open("tags") as fp:
    tags = np.array([l.strip().replace("refs/tags/v", "").split("|") for l in fp.readlines()])
tags_dates = to_datetime(tags[:,1])
order = np.argsort(tags_dates)
skip = {
    "Acts": 101,
    "DD4hep": 61,
}[COMPONENT]
tags = tags[order][skip:]
tags_dates = tags_dates[order][skip:]

with open("eic_container_tags") as fp:
    eic_container_tags = np.array([l.strip().split("|") for l in fp.readlines()])
eic_container_tags_dates = to_datetime(eic_container_tags[:,1])
order = np.argsort(eic_container_tags_dates)
eic_container_tags = eic_container_tags[order]
eic_container_tags_dates = eic_container_tags_dates[order]

with open("result") as fp:
    result = np.array([l.strip().split("|") for l in fp.readlines()])
result_dates = to_datetime(result[:,1])
order = np.argsort(result_dates)
result = result[order]
result_dates = result_dates[order]

assert len(eic_container_tags) == len(result)

result_ixs, tags_ixs = np.nonzero(result[:,0][:,np.newaxis] == tags[:,0][np.newaxis,:])
print(result_ixs)
print(tags_ixs)

def annotate(ax, label, x, y, xytext):
    ax.annotate(label, xy=(x,y),
                xytext=xytext, textcoords='offset points',
                fontsize=5,
                arrowprops={'arrowstyle': '-|>', 'color': 'black', 'linewidth': 1.})

#plt.xkcd()
fig, ax = plt.subplots(figsize=(6.4, 3), layout='constrained')

x = np.arange(len(tags))
ax.plot(tags_dates, x, label=f"{COMPONENT} releases")
ax.plot(eic_container_tags_dates[result_ixs], tags_ixs, label="ePIC stack releases")

ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
ax.xaxis.set_major_locator(mpl.dates.YearLocator())
ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y'))
container_bump_ix = np.nonzero(result[1:,0] != result[:-1,0])[0] + 1
assert np.allclose(result_ixs, range(len(result_ixs)))
for cix, ix in zip(container_bump_ix, tags_ixs[container_bump_ix]):
    if COMPONENT == "Acts":
        annotate(plt.gca(), tags[ix,0], tags_dates[ix], ix, (-200 + ix, ix * 1.2 - 90))
        annotate(plt.gca(), eic_container_tags[cix,0], eic_container_tags_dates[cix], ix, (-100 + ix, -80 + ix * 0.5))
    elif COMPONENT == "DD4hep":
        annotate(plt.gca(), tags[ix,0], tags_dates[ix], ix, (-100 + 2 * ix, ix * 2.4 - 10))
        annotate(plt.gca(), eic_container_tags[cix,0], eic_container_tags_dates[cix], ix, (-10 + 2 * ix, -80 + ix))
plt.ylabel(f"{COMPONENT} version index", loc="top")
plt.legend()
ax.spines[['right', 'top']].set_visible(False)
plt.savefig("fig.png", dpi=300)
plt.show()
