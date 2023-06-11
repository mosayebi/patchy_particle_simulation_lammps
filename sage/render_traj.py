import numpy as np
from pathlib import Path
import os
import subprocess
import argparse


def read_dump(filer, min_step=0, max_step=1e10):
    data = []
    snapshot = {}
    snapshot["traj"] = filer
    cluster_flag = False
    minstep_flag = False
    read_flag = False
    if max_step < min_step:
        min_step, max_step = max_step, min_step
    cnt = 0
    with open(filer, "r") as f:
        while True:
            line = f.readline().strip("\n")
            if not line:
                break
            items = line.split(" ")
            if items[0] == "ITEM:":
                if items[1] == "TIMESTEP":
                    step = int(f.readline().split(" ")[0])
                    if step > max_step:
                        print("max_step reached (%d)" % max_step)
                        print("From %s, last TIMESTEP %d" % (filer, data[-1]["step"]))
                        return data
                    if step >= min_step and minstep_flag == False:
                        print(
                            "From %s, first TIMESTEP reached (%d)" % (filer, min_step)
                        )
                        minstep_flag = True
                    if step >= min_step and step <= max_step:
                        read_flag = True
                        cnt += 1
                        print("%d : reading TIMESTEP %d" % (cnt, step))
                    snapshot["step"] = step
                if items[1] == "NUMBER":
                    N = int(f.readline().split(" ")[0])
                    snapshot["N"] = N
                if items[1] == "BOX":
                    line = f.readline().split(" ")
                    box_x = float(line[1]) - float(line[0])
                    line = f.readline().split(" ")
                    box_y = float(line[1]) - float(line[0])
                    line = f.readline().split(" ")
                    box_z = float(line[1]) - float(line[0])
                    snapshot["box"] = np.array([box_x, box_y, box_z])
                if items[1] == "ATOMS":
                    if len(items) > 7:
                        if "c_cl" in items[7]:
                            cluster_flag = True
                            cluster = np.zeros(N, dtype=np.int)
                        else:
                            cluster = []
                    p_type = np.zeros(N, dtype=np.int)
                    x = np.zeros((N, 3))
                    for i in range(N):
                        line = f.readline().split(" ")
                        p_type[int(line[0]) - 1] = int(line[1])
                        x[int(line[0]) - 1][0] = PBC_wrap(
                            float(line[2]), box_x
                        )  # *float(box_x) - float(box_x)/2
                        x[int(line[0]) - 1][1] = PBC_wrap(
                            float(line[3]), box_y
                        )  # *float(box_y) - float(box_y)/2
                        x[int(line[0]) - 1][2] = PBC_wrap(
                            float(line[4]), box_z
                        )  # *float(box_z) - float(box_z)/2
                        if cluster_flag:
                            cluster[int(line[0]) - 1] = int(line[5])

                    snapshot["p_type"] = p_type
                    snapshot["coords"] = x
                    snapshot["cluster"] = cluster

                    if read_flag:
                        data.append(snapshot.copy())

                    snapshot = {}
                    snapshot["traj"] = filer
    print("From %s, last TIMESTEP %d" % (filer, data[-1]["step"]))
    return data


def conf2tcl(snap, cluster_flag=True, box_flag=True, psi3_flag=False, psi3=[]):
    x = snap["coords"]
    p_type = snap["p_type"]
    box = snap["box"]
    step = snap["step"]
    N = snap["N"]
    cluster = snap["cluster"]
    if cluster_flag and len(cluster) == 0:
        print("Warning: no cluster data found. Cluster coloring is off.")
        cluster_flag = False
        psi3_flag = False
    if len(psi3) != 0:
        psi3_flag = True
        cluster_flag = False
        psi3 = np.nan_to_num(psi3)

    ret = "color Display Background white\n"
    ret += "display projection orthographic\n"
    ret += "axes location off\n"
    ret += "mol new\n"
    ret += "color scale method RGB\n"
    ret += "graphics 0 delete all\n"

    if box_flag:
        box_radius = 0.1
        ret += "graphics 0 color 0\n"
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                -box[0] / 2.0,
                -box[1] / 2.0,
                -box[2] / 2.0,
                box[0] / 2.0,
                -box[1] / 2.0,
                -box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                -box[0] / 2.0,
                -box[1] / 2.0,
                box[2] / 2.0,
                box[0] / 2.0,
                -box[1] / 2.0,
                box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                -box[0] / 2.0,
                +box[1] / 2.0,
                -box[2] / 2.0,
                -box[0] / 2.0,
                -box[1] / 2.0,
                -box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                -box[0] / 2.0,
                +box[1] / 2.0,
                box[2] / 2.0,
                -box[0] / 2.0,
                -box[1] / 2.0,
                box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                box[0] / 2.0,
                +box[1] / 2.0,
                box[2] / 2.0,
                box[0] / 2.0,
                -box[1] / 2.0,
                box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                box[0] / 2.0,
                +box[1] / 2.0,
                -box[2] / 2.0,
                box[0] / 2.0,
                -box[1] / 2.0,
                -box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                -box[0] / 2.0,
                +box[1] / 2.0,
                -box[2] / 2.0,
                box[0] / 2.0,
                +box[1] / 2.0,
                -box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                -box[0] / 2.0,
                +box[1] / 2.0,
                box[2] / 2.0,
                box[0] / 2.0,
                +box[1] / 2.0,
                box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                -box[0] / 2.0,
                -box[1] / 2.0,
                -box[2] / 2.0,
                -box[0] / 2.0,
                -box[1] / 2.0,
                +box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                box[0] / 2.0,
                -box[1] / 2.0,
                -box[2] / 2.0,
                box[0] / 2.0,
                -box[1] / 2.0,
                +box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                box[0] / 2.0,
                box[1] / 2.0,
                -box[2] / 2.0,
                box[0] / 2.0,
                box[1] / 2.0,
                +box[2] / 2.0,
            )
        )
        ret += (
            "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 20 filled yes\n"
            % (
                -box[0] / 2.0,
                box[1] / 2.0,
                -box[2] / 2.0,
                -box[0] / 2.0,
                box[1] / 2.0,
                +box[2] / 2.0,
            )
        )

    if psi3_flag:
        cmax = 1.0
        cmin = 0.0
        # if len(psi3[np.nonzero(np.nan_to_num(psi3))]) > 0 :
        #    cmin = np.min(psi3[np.nonzero(psi3)])
        #    if cmin == 1 :
        #     cmin = 0
        # print cmin, cmax

    # N=96
    ret += "graphics 0 color 7\n"  # blue
    for i in range(N):
        if i % 16 == 0 and p_type[i] == 1:
            if cluster_flag:
                ret += vmd_hub(
                    x[i, :],
                    x[i + 2, :],
                    box,
                    int(float(cluster[i]) / max(cluster[:]) * 1023 + 1),
                )
            elif psi3_flag:
                c = np.nan_to_num(psi3[get_mol_id(i)]) ** 2
                if c == 0:
                    c = 1
                else:
                    c = int((c - cmin) / (cmax - cmin) * 950 + 70)
                ret += vmd_hub(x[i, :], x[i + 2, :], box, c)
            else:
                ret += vmd_hub(x[i, :], x[i + 2, :], box, 0)
        else:
            continue

    ret += "graphics 0 color 0\n"  # blue
    for i in range(N):
        # ret += '#%s  %s %s %s %s\n'% (i,p_type[i], x[i,0],x[i,1],x[i,2])
        if i % 16 == 10 and p_type[i + 3] == 8:
            if cluster_flag:
                ret += vmd_hub(
                    x[i, :],
                    x[i + 2, :],
                    box,
                    int(float(cluster[i]) / max(cluster[:]) * 1023 + 1),
                )
            elif psi3_flag:
                c = np.nan_to_num(psi3[get_mol_id(i)]) ** 2
                if c == 0:
                    c = 1
                else:
                    c = int((c - cmin) / (cmax - cmin) * 950 + 70)
                ret += vmd_hub(x[i, :], x[i + 2, :], box, c)
            else:
                ret += vmd_hub(x[i, :], x[i + 2, :], box, 0)
        else:
            continue

    ret += "graphics 0 color 1\n"  # red
    for i in range(N):
        if i % 16 == 10 and p_type[i + 3] == 11:
            if cluster_flag:
                ret += vmd_hub(
                    x[i, :],
                    x[i + 2, :],
                    box,
                    int(float(cluster[i]) / max(cluster[:]) * 1023 + 1),
                )
            elif psi3_flag:
                c = np.nan_to_num(psi3[get_mol_id(i)]) ** 2
                if c == 0:
                    c = 1
                else:
                    c = int((c - cmin) / (cmax - cmin) * 950 + 70)
                ret += vmd_hub(x[i, :], x[i + 2, :], box, c)
            else:
                ret += vmd_hub(x[i, :], x[i + 2, :], box, 0)
        else:
            continue

    ret += "graphics 0 color 6\n"  # silver
    for i in range(N):
        if i % 16 == 10 and p_type[i] == 1:
            if cluster_flag:
                ret += vmd_bond(
                    x[i, :],
                    x[i - 10, :],
                    box,
                    int(float(cluster[i]) / max(cluster[:]) * 1023 + 1),
                )
            elif psi3_flag:
                c = np.nan_to_num(psi3[get_mol_id(i - 10)]) ** 2
                if c == 0:
                    c = 1
                else:
                    c = int((c - cmin) / (cmax - cmin) * 950 + 70)
                ret += vmd_bond(x[i, :], x[i - 10, :], box, c)
            else:
                ret += vmd_bond(x[i, :], x[i - 10, :], box, 0)
        else:
            continue

    ret += "display resetview\n"
    ret += "rotate x by 10\n"
    ret += "rotate y by -10\n"

    return ret


def vmd_bond(x1, x2, box, c=0):
    n = PBC(x2 - x1, box)
    r1 = x1
    r2 = r1 + n
    ret = ""
    if not c == 0:
        ret += "graphics 0 color %s\n" % c
    ret += (
        "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.1 resolution 30 filled yes\n"
        % (r1[0], r1[1], r1[2], r2[0], r2[1], r2[2])
    )
    return ret


def vmd_hub(x1, x2, box, c=0):
    n = PBC(x2 - x1, box)
    mag_n = np.linalg.norm(n)
    r1 = x1 - 0.5 * n / mag_n
    r2 = r1 + 3 * n / mag_n
    ret = ""
    if not c == 0:
        ret += "graphics 0 color %s\n" % c
    ret += (
        "graphics 0 cylinder {%lf %lf %lf} {%lf %lf %lf} radius 0.4 resolution 30 filled yes\n"
        % (r1[0], r1[1], r1[2], r2[0], r2[1], r2[2])
    )
    return ret


def vmd_patch(x, c):
    ret = ""
    if not c == 0:
        ret += "graphics 0 color %s\n" % c
    ret += "graphics 0 sphere {%lf %lf %lf} radius 0.2 resolution 30\n" % (
        x[0],
        x[1],
        x[2],
    )
    return ret


def write_tcl(tcl_out, filename="out.tcl"):
    try:
        f = open(filename, "w")
    except:
        print("tcl_output: cannot open output file. Dying")
    f.write(tcl_out)


def render_tcl_file(
    res_x,
    res_y,
    tcl_out_file="out.tcl",
    png_file="out.png",
    tachyon_path=Path(
        "/Applications/VMD 1.9.4a57-x86_64-Rev12.app/Contents/vmd/tachyon_MACOSXX86_64"
    ),
    vmd_path=Path(
        "/Applications/VMD 1.9.4a57-x86_64-Rev12.app/Contents/MacOS/startup.command"
    ),
):
    assert tachyon_path.exists(), "tachyon was not found in `{tachyon_path}`"
    assert vmd_path.exists(), "vmd was not found in `{tachyon_path}`"

    tcl_out = "source %s\nrender Tachyon out.render\n" % (tcl_out_file)
    vmdin = os.popen(f'"{vmd_path}" -dispdev none', "w")
    vmdin.write("%s" % tcl_out)
    vmdin.flush()

    tga_file = str(png_file) + ".tga"
    args = "-aasamples 24 out.render -format TGA -res %s %s -o %s" % (
        res_x,
        res_y,
        tga_file,
    )
    command = []
    command.append(str(tachyon_path))
    command = command + args.split()
    #print(command)
    exe = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = exe.communicate()
    #print(stdout, stderr)


def PBC_wrap(x, box):
    while x < -box * 0.5:
        x = x + box
    while x >= box * 0.5:
        x = x - box
    return x


def PBC(d, box):
    d[0] = d[0] - np.rint(d[0] / box[0]) * box[0]
    d[1] = d[1] - np.rint(d[1] / box[1]) * box[1]
    d[2] = d[2] - np.rint(d[2] / box[2]) * box[2]
    return d


def run(args):
    traj_data = read_dump(args.lammps_dump, args.min_timestep, args.max_timestep)
    traj_tcl = [
        conf2tcl(
            snap,
            cluster_flag=args.cluster_flag,
            box_flag=args.box_flag,
            psi3_flag=False,
            psi3=[],
        )
        for s, snap in enumerate(traj_data)
    ]

    args.output.mkdir(exist_ok=True, parents=True)
    for snap, tcl in zip(traj_data, traj_tcl):
        tcl_file = args.output / f"{args.lammps_dump.name}_{snap['step']}.tcl"
        write_tcl(tcl, filename=tcl_file)
        if args.render:
            print(["rendering tcl files ..."])
            png_file = args.output / f"{args.lammps_dump.name}_{snap['step']}.png"
            render_tcl_file(1024, 1024, tcl_out_file=tcl_file, png_file=png_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("lammps_dump", type=Path)
    parser.add_argument("--output", type=Path, default=Path("render_output"))
    parser.add_argument("--min_timestep", default=0, type=int)
    parser.add_argument("--max_timestep", default=1e9, type=int)
    parser.add_argument("--render", action="store_true", default=False)
    parser.add_argument("--cluster_flag", action="store_true", default=False)
    parser.add_argument("--box_flag", action="store_true", default=True)
    parser.add_argument(
        "--tachyon_path",
        type=Path,
        default=Path(
            "/Applications/VMD 1.9.4a57-x86_64-Rev12.app/Contents/vmd/tachyon_MACOSXX86_64"
        ),
    )
    parser.add_argument(
        "--vmd_path",
        type=Path,
        default=Path(
            "/Applications/VMD 1.9.4a57-x86_64-Rev12.app/Contents/MacOS/startup.command"
        ),
    )
    args = parser.parse_args()

    # check for vmd paths, if args.render
    if args.render:
        print(f"[For rendering, vmd/tachyon is required.]")
        print(
            f"default paths are:\n {dict(vmd_path=args.vmd_path, tachyon_path=args.tachyon_path)}"
        )
        assert (
            args.tachyon_path.exists()
        ), "tachyon was not found in `{args.tachyon_path}`. Please pass the correct path via --tachyon_path"
        assert (
            args.vmd_path.exists()
        ), "vmd was not found in `{tachyon_path}`. Please pass the correct path via --vmd_path"

    run(args)
