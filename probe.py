import argparse


import MDAnalysis as mda
from ggmolvis.ggmolvis import GGMolVis
import bpy
import yaml


def set_background(rgba: tuple):
    # NOTE: manual bpy usage--this functionality should be upstreamed to ggmolvis
    bpy.context.scene.render.film_transparent = False
    world = bpy.context.scene.world
    if not world.use_nodes:
        world.use_nodes = True
    node_tree = world.node_tree
    bg_node = node_tree.nodes.get("Background") or node_tree.nodes.new("ShaderNodeBackground")
    bg_node.inputs["Color"].default_value = rgba
    output_node = node_tree.nodes.get("World Output") or node_tree.nodes.new("ShaderNodeOutputWorld")
    if not any(link.to_node == output_node for link in bg_node.outputs[0].links):
        node_tree.links.new(bg_node.outputs[0], output_node.inputs[0])


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
    mode = p_config["mode"]
    background_color = p_config["background_color"]
    focal_length = p_config["focal_length"]

    if mode == "movie":
        outfile_suffix = "mp4"
    else:
        outfile_suffix = "png"


    if blender_engine == "CYCLES":
        # NOTE: manual usage of bpy to set blender device, should upstream to ggmolvis
        bpy.context.scene.cycles.device = blender_engine_device

    ggmv = GGMolVis()
    u = mda.Universe(topology_path, trajectory_path)
    system = u.select_atoms(sel_string)
    all_mol = ggmv.molecule(system, lens=focal_length)

    # NOTE: manual bpy usage here for setting Blender engine--should be abstracted to ggmolvis
    bpy.context.scene.render.engine = blender_engine

    # NOTE: manual bpy usage--need merge of frame handling capability from i.e.,
    # https://github.com/yuxuanzhuang/ggmolvis/pull/16
    bpy.context.scene.frame_start = frame_start
    bpy.context.scene.frame_end = frame_end

    set_background(rgba=(background_color["red"],
                         background_color["green"],
                         background_color["blue"],
                         background_color["alpha"]))

    all_mol.render(resolution=(blender_resolution_x, blender_resolution_y),
                   filepath=f"{render_filename}.{outfile_suffix}",
                   mode=mode)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="North Star rendering suite for probing ggmolvis capabilities/gaps")
    parser.add_argument("render_config_file", help="Path to the render configuration YML file.")
    args = parser.parse_args()
    with open(args.render_config_file) as config_file:
        p_config = yaml.load(config_file, Loader=yaml.FullLoader)
    main(p_config=p_config)
