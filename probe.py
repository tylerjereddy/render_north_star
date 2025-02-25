import argparse


import MDAnalysis as mda
from ggmolvis.ggmolvis import GGMolVis
import bpy
import yaml


def main(p_config: dict):
    """
    Parameters
    ----------
    p_config: dict
        Dictionary of rendering parameters parsed from the YAML
        input deck.
    """
    # read in the YAML input deck parameters
    topology_path = p_config["topology_path"]
    trajectory_path = p_config["trajectory_path"]
    blender_engine = p_config["blender_engine"]
    blender_engine_device = p_config["blender_engine_device"]
    blender_resolution_x = p_config["blender_resolution_x"]
    blender_resolution_y = p_config["blender_resolution_y"]
    render_filename = p_config["render_filename"]
    frame_start = p_config["frame_start"]
    frame_end = p_config["frame_end"]
    sel_string = p_config["sel_string"]


    if blender_engine == "CYCLES":
        # NOTE: manual usage of bpy to set blender device, should upstream to ggmolvis
        bpy.context.scene.cycles.device = blender_engine_device

    ggmv = GGMolVis()
    u = mda.Universe(topology_path, trajectory_path)
    system = u.select_atoms(sel_string)
    all_mol = ggmv.molecule(system)

    # NOTE: manual bpy usage here for setting Blender engine--should be abstracted to ggmolvis
    bpy.context.scene.render.engine = blender_engine

    # NOTE: manual bpy usage--need merge of frame handling capability from i.e.,
    # https://github.com/yuxuanzhuang/ggmolvis/pull/16
    bpy.context.scene.frame_start = frame_start
    bpy.context.scene.frame_end = frame_end

    all_mol.render(resolution=(blender_resolution_x, blender_resolution_y),
                   filepath=render_filename,
                   mode="movie")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="North Star rendering suite for probing ggmolvis capabilities/gaps")
    parser.add_argument("render_config_file", help="Path to the render configuration YML file.")
    args = parser.parse_args()
    with open(args.render_config_file) as config_file:
        p_config = yaml.load(config_file, Loader=yaml.FullLoader)
    main(p_config=p_config)
