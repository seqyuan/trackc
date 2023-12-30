import sys

import click
import pandas as pd
import yaml

import trackc as tc


def _get_yaml_data(yaml_file):
    file = open(yaml_file, "r", encoding="utf-8")
    file_data = file.read()
    file.close()

    data = yaml.load(file_data)
    return data


def init_ax(ax_dic, height=1, hspace=0.06):
    if "height" in ax_dic.keys():
        if ax_dic["height"] != None:
            height = float(ax_dic["height"])

    if "hspace" in ax_dic.keys():
        if ax_dic["hspace"] != None:
            hspace = float(ax_dic["hspace"])

    return height, hspace


def _make_spec(conf):
    df = pd.DataFrame({"ax": [tt["ax"] for i, tt in enumerate(conf["trackc"])]})
    df["track_type"] = [tt["track_type"] for i, tt in enumerate(conf["trackc"])]
    axs = pd.DataFrame(
        {"rank_u": range(len(df["ax"].unique()))}, index=df["ax"].unique()
    )
    axs["height"] = 1
    axs["hspace"] = 0.05
    axs["bottom"] = "bottom"

    for i, tt in enumerate(conf["trackc"]):
        height, hspace = init_ax(
            tt, axs.loc[tt["ax"], "height"], axs.loc[tt["ax"], "hspace"]
        )
        axs.loc[tt["ax"], "height"] = height
        axs.loc[tt["ax"], "hspace"] = hspace

    return df, axs


def hicContactMap(Paras, ax):
    mapc_paras = Paras["mapC"]
    mat = None
    mat2 = None
    if "mat" in Paras.keys():
        if Paras["mat"]["method"] == "extractCisContact":
            del Paras["mat"]["method"]
            mat = tc.tl.extractCisContact(**Paras["mat"])

        elif Paras["mat"]["method"] == "extractContactRegions":
            if (
                "row_regions" in Paras["mat"].keys()
                or "col_regions" in Paras["mat"].keys()
            ):
                del Paras["mat"]["method"]
                mat_obj = tc.tl.extractContactRegions(
                    clr=Paras["mat"]["clr"], row_regions=region
                )
                mat = mat_obj.cmat
        else:
            pass

    if "mat2" in Paras.keys():
        if Paras["mat2"]["method"] == "extractCisContact":
            del Paras["mat2"]["method"]
            mat = tc.tl.extractCisContact(**Paras["mat2"])

        elif Paras["mat2"]["method"] == "extractContactRegions":
            if (
                "row_regions" in Paras["mat2"].keys()
                or "col_regions" in Paras["mat2"].keys()
            ):
                del Paras["mat2"]["method"]
                mat_obj = tc.tl.extractContactRegions(
                    clr=Paras["mat2"]["clr"], row_regions=region
                )
                mat = mat_obj.cmat
        else:
            pass

    if mat is None:
        tc.pl.mapC(ax=ax, mat2=mat2, **Paras["mapC"])
    else:
        if mat2 is not None:
            tc.pl.mapC(ax=ax, mat=mat, mat2=mat2, **Paras["mapC"])
        else:
            tc.pl.mapC(ax=ax, mat=mat, **Paras["mapC"])


@click.command()
@click.argument("config", metavar="<trackc-conf.yml>")
@click.option(
    "--regions",
    "-r",
    default=None,
    help="genome regions, \
              eg.:'18:47950000-48280000 18:75280000-74850000'",
)
@click.option(
    "--outfile", "-o", default="trackc.pdf", help="output filename. default=trackc.pdf"
)
@click.option(
    "--basefigsize",
    "-s",
    default="6,1",
    help='base figsize: width,height. default="6,1"\
              The height option in the config.yml is relative to the base figsize height',
)
def cli(config, regions, outfile, basefigsize):
    conf = _get_yaml_data(config)
    fs = basefigsize.split(",")
    bfs = (float(fs[0]), float(fs[1]))
    if regions != None:
        region_tmp = regions.split(" ")
        if len(region_tmp) > 1:
            regions = region_tmp

    # make grid spec
    df, axs_pre = _make_spec(conf)

    ten = tc.tenon(figsize=bfs)
    for i, row in axs_pre.iterrows():
        ten.add(pos=row["bottom"], height=row["height"], hspace=row["hspace"])

    regions_type = [
        "multi_scale_track",
        "bed_track",
        "bw_track",
        "virtual4C",
        "links_track",
        "gene_track",
    ]
    region_type = ["scale_track"]

    for i, tt in enumerate(conf["trackc"]):
        track_type = tt["track_type"]
        if track_type == None:
            ax = ten.axs(axs_pre.loc[tt["ax"], "rank_u"])
            ax.set_axis_off()
            continue
        paras = tt["track_para"]
        if track_type in regions_type:
            if "regions" not in paras.keys():
                if regions != None:
                    paras["regions"] = regions
                else:
                    print("region not set")
                    sys.exit(1)

        if track_type in region_type:
            if "region" not in paras.keys():
                if regions != None:
                    paras["region"] = regions
                else:
                    # print('region para not set in cli and conf.yml')
                    pass

        if track_type == "hicmap":
            if "mat" in paras.keys():
                if paras["mat"]["method"] == "extractCisContact":
                    if "region" not in paras["mat"].keys():
                        paras["mat"]["region"] = regions
                elif paras["mat"]["method"] == "extractContactRegions":
                    if (
                        "row_regions" not in paras["mat"].keys()
                        and "col_regions" not in paras["mat"].keys()
                    ):
                        paras["mat"]["row_regions"] = regions
                else:
                    pass

            if "mat2" in paras.keys():
                if paras["mat2"]["method"] == "extractCisContact":
                    if "region" not in paras["mat"].keys():
                        paras["mat2"]["region"] = regions
                elif paras["mat2"]["method"] == "extractContactRegions":
                    if (
                        "row_regions" not in paras["mat2"].keys()
                        and "col_regions" not in paras["mat2"].keys()
                    ):
                        paras["mat2"]["row_regions"] = regions
                else:
                    pass

        ax = ten.axs(axs_pre.loc[tt["ax"], "rank_u"])

        if track_type == "multi_scale_track":
            tc.pl.multi_scale_track(**paras, ax=ax)

        if track_type == "scale_track":
            tc.pl.scale_track(**paras, ax=ax)

        if track_type == "bed_track":
            tc.pl.bed_track(**paras, ax=ax)

        if track_type == "bw_track":
            tc.pl.bw_track(**paras, ax=ax)

        if track_type == "virtual4C":
            tc.pl.virtual4C(**paras, ax=ax)

        if track_type == "links_track":
            tc.pl.links_track(**paras, ax=ax)

        if track_type == "gene_track":
            tc.pl.gene_track(**paras, ax=ax)

        if track_type == "hicmap":
            hicContactMap(paras, ax=ax)

        if track_type == "zoomin":
            tc.pl.zoomin(**paras, ax=ax)

    tc.savefig(outfile)
