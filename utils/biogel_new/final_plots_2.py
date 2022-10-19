from turtle import Turtle
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import matplotlib as mpl
import cmasher as cmr
import pandas as pd
from matplotlib import font_manager
from matplotlib import gridspec

fe = font_manager.FontEntry(
    fname="/Users/ada/Library/Fonts/Lato-Light.ttf", name="Lato Light"
)
font_manager.fontManager.ttflist.insert(0, fe)
mpl.rcParams["font.family"] = fe.name

plt.style.use("./publication.mplstyle")


def parse_file(file):
    df = pd.read_pickle(file)
    df_param = pd.read_csv("parameters.txt", delim_whitespace=True, dtype="str")
    df_param.rename(columns={"ID": "state_id"}, inplace=True)
    df_param = df_param.astype(
        {
            "N": float,
            "rho": float,
            "state_id": str,
            "kb": float,
            "plink": float,
            "epsilon": float,
        }
    )
    Nmax = 51000

    df_param["Nparticles"] = (Nmax * df_param.rho) / df_param.N

    df = pd.merge(df, df_param, on="state_id")
    df = df.astype(
        {
            "state_id": str,
            "time": int,
            "N": float,
            "kb": float,
            "plink": float,
            "rho": float,
            "epsilon": float,
            "largest_domain": float,
            "node_connectivity": int,
            "pstart": str,
        }
    )

    return df, df_param


def filter_conditions(df, x, conditions):
    for key, val in conditions.items():
        if key != x:
            df = df[df[key] == val]

    return df


def evaluate_distri_degree(df, df_param):
    dg = pd.DataFrame(
        {"degree_sequence": df.groupby("state_id").degree_sequence.apply(list)}
    ).reset_index()
    dg.degree_sequence = dg.degree_sequence.apply(
        lambda x: [item for sublist in x for item in sublist]
    )
    dg = pd.merge(dg, df_param, on="state_id")

    return dg


def evaluate_distri_domain_size(df, df_param):
    dg = pd.DataFrame(
        {"domain_lengths": df.groupby("state_id").domain_lengths.apply(list)}
    ).reset_index()

    dg.domain_lengths = dg.domain_lengths.apply(
        lambda x: [item for sublist in x for item in sublist]
    )
    dg = pd.merge(dg, df_param, on="state_id")

    return dg


def evaluate_mean_largest_cluster(df, df_param):
    df["percentage_largest"] = df.largest_domain / df.Nparticles
    dg = pd.DataFrame(
        {
            "mean_percentage_largest": df.groupby("state_id").percentage_largest.mean(),
            "std_percentage_largest": df.groupby("state_id").percentage_largest.std(),
        }
    ).reset_index()
    dg = pd.merge(dg, df_param, on="state_id")
    return dg


def plot_distri_domain_size(df, p_params, ax=None):

    if ax is None:
        ax = plt.gca()

    x = p_params["var"]
    if p_params["log_scale"] == True:
        ax.set_yscale("log")

    for lx, x_i in enumerate(sorted(df[x].unique())):
        dg = df[(df[x] == x_i)]
        arr = dg["domain_lengths"].values
        bins = int(np.floor(np.max(arr[0]) / p_params["bin_scale"]))

        if p_params["normed"] == True:
            arr[0] = arr[0] / dg["Nparticles"].values

        hist, bin_edges = np.histogram(
            np.array(arr[0]),
            bins=bins,
            density=True,
        )
        xarr = bin_edges[:-1]

        label = "$\{}$".format(x)
        ax.plot(
            xarr,
            hist,
            marker=None,
            ms=0,
            label="{} = {}".format(label, x_i),
        )

    ax.set_xlabel(p_params["xlabel"])
    ax.set_ylabel(p_params["ylabel"])
    ax.set_xlim((p_params["xmin"], p_params["xmax"]))
    ax.set_ylim((p_params["ymin"], p_params["ymax"]))
    plt.tight_layout()


def plot_distri_degree(df, plot_params, ax=None):

    if ax is None:
        ax = plt.gca()

    x = p_params["var"]
    for lx, x_i in enumerate(sorted(df[x].unique())):
        dg = df[(df[x] == x_i)]
        arr = dg["degree_sequence"].values
        bins = np.max(arr[0])

        hist, bin_edges = np.histogram(np.array(arr[0]), bins=bins, density=True)
        xarr = bin_edges[:-1]

        label = "$\{}$".format(x)
        ax.plot(
            xarr,
            hist,
            label="{} = {}".format(label, x_i),
        )
        ax.set_xlabel(p_params["xlabel"])
        ax.set_ylabel(p_params["ylabel"])
        ax.set_xlim((p_params["xmin"], p_params["xmax"]))
        ax.set_ylim((p_params["ymin"], p_params["ymax"]))
        plt.tight_layout()


def plot_mean_largest_cluster(df, p_params, ax=None):

    if ax is None:
        ax = plt.gca()

    x = p_params["var"]
    arr = df[["mean_percentage_largest", "std_percentage_largest", x]].values
    arr = arr[arr[:, 2].argsort()]

    ax.errorbar(arr[:, 2], arr[:, 0], arr[:, 1], c="k", alpha=0.5)

    label = "$\{}$".format(x)
    for i, x_i in enumerate(sorted(df[x].unique())):
        ax.plot(
            arr[i, 2],
            arr[i, 0],
            label="{} = {}".format(label, x_i),
            zorder=20,
        )

    # add line for expected number of polymers without crosslinker
    # P = exp(-p) p = Nbeads*plink
    Nbeads = 30
    plink = 0.1
    x = np.arange(1, 6)
    y = np.ones(5) * (1 - np.exp(-Nbeads * plink))
    ax.plot(x, y, color="r", linestyle="--")
    ax.set_xlabel(p_params["xlabel"])
    ax.set_ylabel(p_params["ylabel"])
    ax.set_ylim((p_params["ymin"], p_params["ymax"]))
    ax.set_xlim((p_params["xmin"], p_params["xmax"]))


def plot_gr(df, ax):

    if ax == None:
        ax = plt.gca()

    arr = df.values
    epsi = [1, 1.5, 2, 3, 5]
    ax.set_yscale("log")

    for i in range(5):
        ax.plot(
            arr[:, i * 2],
            arr[:, i * 2 + 1],
            label="{} = {}".format("$\epsilon$", epsi[i]),
            ms=0,
        )

    ax.set_xlabel("$r /\sigma$")
    ax.set_ylabel("$g(r /\sigma$)")
    ax.set_ylim((0.3, 100))
    ax.set_xlim((0.7, 5))


def plot_nbonded(df, ax):
    arr = df[["epsilon", "nlinked"]].values
    ax.plot(arr[:, 0], arr[:, 1], c="r")
    ax.set_xlim((0.9, 5.1))
    ax.set_ylim((0, 1))
    ax.set_xlabel("$\epsilon$")
    ax.set_ylabel("$n_{l}$")


def plot_mean_pe(df, conditions, ax):

    if ax == None:
        ax = plt.gca()

    n = len(df["rho"].unique()) + 1
    clist = [plt.get_cmap("cmr.fall_r")(1.0 * i / n) for i in range(n)]

    for i, rho_x in enumerate(sorted(df["rho"].unique())):
        conditions["rho"] = rho_x
        df_filtered = filter_conditions(df, "epsilon", conditions)
        arr = df_filtered[["epsilon", "mean_pe"]].values
        ax.plot(
            arr[:, 0], arr[:, 1], c=clist[i + 1], label="$\\rho = {}$".format(rho_x)
        )

    ax.set_xlabel("$\epsilon$")
    ax.set_ylabel("$E_{l} / [ \epsilon \cdot n_{b} ]$")
    ax.legend(ncol=1, bbox_to_anchor=(1, 1), loc="upper right")
    plt.tight_layout()


def plot_pe(df, ax):
    if ax == None:
        ax = plt.gca()

    arr = df.values
    epsilons = [1, 1.5, 2, 3, 5]
    ax.set_yscale("log")
    for i, ei in enumerate(epsilons):
        ax.plot(
            arr[:, i * 2], arr[:, i * 2 + 1], ms=0, label="$\epsilon = {}$".format(ei)
        )
    ax.set_xlim((-1.0, 0))
    ax.set_ylim((0.1, 50))
    ax.set_xlabel("$E_{l}/\epsilon$")
    ax.set_ylabel("P")
    plt.tight_layout()


def plot_local_density(df_dict, plot_params, ax):

    if ax == None:
        ax = plt.gca()

    ax.set_yscale("log")

    for epsi in df_dict:
        arr = df_dict[epsi].values
        ax.plot(arr[:, 0], arr[:, 1], ms=0, label="$\epsilon = {}$".format(epsi))

    ax.set_xlim((plot_params["xmin"], plot_params["xmax"]))
    ax.set_ylim((plot_params["ymin"], plot_params["ymax"]))
    ax.set_xlabel("$\\rho_{l}$")
    ax.set_ylabel("P")
    plt.tight_layout()


if __name__ == "__main__":
    n = 5
    color_list = [plt.get_cmap("cmr.lavender_r")(1.0 * i / n) for i in range(n)]
    default_cycler = cycler(color=color_list)
    plt.rc("axes", prop_cycle=default_cycler)

    # Fig 2a,b,c - polymer properties std conditions
    df, df_param = parse_file("network_data_polymer_all_runs.pickle")
    std_conditions = {
        "epsilon": 5,
        "rho": 0.2,
        "kb": 30,
        "N": 30,
        "plink": 0.1,
    }
    df_filtered = filter_conditions(df, "epsilon", std_conditions)

    df_domain_size = evaluate_distri_domain_size(df_filtered, df_param)
    df_polymer_degree = evaluate_distri_degree(df_filtered, df_param)
    df_mean_largest = evaluate_mean_largest_cluster(df_filtered, df_param)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    # degree plot
    p_params = {
        "xmin": 0,
        "xmax": 30,
        "ymin": 0,
        "ymax": 0.55,
        "xlabel": "polymer degree",
        "ylabel": "P",
        "var": "epsilon",
    }

    plot_distri_degree(df_polymer_degree, p_params, ax1)

    # domain size
    p_params = {
        "xmin": 0,
        "xmax": 1,
        "ymin": 0.01,
        "ymax": 35,
        "xlabel": "fractional polymer domain size",
        "ylabel": "P",
        "var": "epsilon",
        "bin_scale": 10,
        "log_scale": True,
        "normed": True,
    }

    plot_distri_domain_size(df_domain_size, p_params, ax2)

    # percolation
    p_params = {
        "xmin": 0.9,
        "xmax": 5.1,
        "ymin": 0,
        "ymax": 1,
        "xlabel": "$\epsilon$",
        "ylabel": "fraction of largest polymer cluster",
        "var": "epsilon",
    }

    plot_mean_largest_cluster(df_mean_largest, p_params, ax3)

    plt.legend(ncol=1, bbox_to_anchor=(1.02, 0.5), loc="upper right", fontsize=16)
    plt.savefig("fig_2abc.pdf")

    # Fig 3a,b,c - crosslinker properties for std conditions
    df, df_param = parse_file("network_data_crosslinker_all_runs.pickle")
    std_conditions = {
        "epsilon": 5,
        "rho": 0.2,
        "kb": 30,
        "N": 30,
        "plink": 0.1,
    }
    df_filtered = filter_conditions(df, "epsilon", std_conditions)
    df_crosslinker_degree = evaluate_distri_degree(df_filtered, df_param)
    df_crosslinker_domain_size = evaluate_distri_domain_size(df_filtered, df_param)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    # degree plot
    p_params = {
        "xmin": 0,
        "xmax": 15,
        "ymin": 0,
        "ymax": 1,
        "xlabel": "crosslinker degree",
        "ylabel": "P",
        "var": "epsilon",
    }

    # cluster size plot
    plot_distri_degree(df_crosslinker_degree, p_params, ax1)

    p_params = {
        "xmin": 0,
        "xmax": 60,
        "ymin": 0,
        "ymax": 0.3,
        "xlabel": "crosslinker cluster size",
        "ylabel": "P",
        "var": "epsilon",
        "bin_scale": 1,
        "log_scale": False,
        "normed": False,
    }
    plot_distri_domain_size(df_crosslinker_domain_size, p_params, ax2)

    # g(r) plot
    df_gr = pd.read_csv("Fig2-grPs3.dat", header=None, delim_whitespace=True)
    plot_gr(df_gr, ax3)

    # inset to ax1 number of bonded linkers

    plt.legend(ncol=1, bbox_to_anchor=(1.02, 0.75), loc="upper right", fontsize=16)

    df = pd.read_csv(
        "Fig3-Nlinked.dat",
        header=None,
        names=["epsilon", "rho", "kb", "N", "plink", "nlinked"],
        delim_whitespace=True,
    )

    df_filtered = filter_conditions(df, "epsilon", std_conditions)
    axins = fig.add_axes([0.13, 0.5, 0.15, 0.3])
    plot_nbonded(df_filtered, axins)

    plt.savefig("fig_3abc.pdf")

# Fig 4abc mean potential and potential energy distribtuions per link and particle

# fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

fig = plt.figure(figsize=(10, 5))

ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 4)

# mean potential
df = pd.read_csv(
    "Fig5a-epot2.dat",
    header=None,
    names=["epsilon", "rho", "kb", "N", "plink", "mean_pe"],
    delim_whitespace=True,
)

plot_mean_pe(df, std_conditions, ax1)

# ax1.legend(ncol=1, bbox_to_anchor=(1.02, 0.75), loc="upper right", fontsize=16)

# mean potential probability for rho = 0.05

df = pd.read_csv(
    "Fig5c-pepotB8Ps12.dat",
    header=None,
    names=["xe1", "pe1", "xe2", "pe2", "xe3", "pe3", "xe4", "pe4", "xe5", "pe5"],
    delim_whitespace=True,
)
plot_pe(df, ax2)

df = pd.read_csv(
    "Fig5b-pepotB8Ps3.dat",
    header=None,
    names=["xe1", "pe1", "xe2", "pe2", "xe3", "pe3", "xe4", "pe4", "xe5", "pe5"],
    delim_whitespace=True,
)

plot_pe(df, ax3)
plt.legend(ncol=2, bbox_to_anchor=(0.8, 0.9), loc="upper right")
plt.savefig("fig_4abc.pdf")

# Fig local density

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

epsilon = [1, 1.5, 2, 3, 5]
data_files = [
    "Local2sigma.1.dat",
    "Local2sigma.2.dat",
    "Local2sigma.3.dat",
    "Local2sigma.4.dat",
    "Local2sigma.5.dat",
]

df_dict = {}
for i, epsi in enumerate(epsilon):
    df = pd.read_csv(
        data_files[i],
        header=None,
        names=["x", "p_local"],
        delim_whitespace=True,
    )
    df_dict[epsi] = df

plot_params = {
    "ymin": 1e-8,
    "ymax": 1e-1,
    "xmin": 0,
    "xmax": 1.2,
}

plot_local_density(df_dict, plot_params, ax1)

data_files = [
    "Local15sigma.1.dat",
    "Local15sigma.2.dat",
    "Local15sigma.3.dat",
    "Local15sigma.4.dat",
    "Local15sigma.5.dat",
]

df_dict = {}
for i, epsi in enumerate(epsilon):
    df = pd.read_csv(
        data_files[i],
        header=None,
        names=["x", "p_local"],
        delim_whitespace=True,
    )
    df_dict[epsi] = df

plot_params = {
    "ymin": 1e-8,
    "ymax": 5e-5,
    "xmin": 0.15,
    "xmax": 0.23,
}

plot_local_density(df_dict, plot_params, ax2)
plt.legend(ncol=1, bbox_to_anchor=(0.4, 0.9), loc="upper right")
plt.savefig("fig_5ab.pdf")
